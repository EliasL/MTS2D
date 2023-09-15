#include <mkl.h>

#define ae_bool  char
#define ae_true  1
#define ae_false 0
#define ae_int_t ptrdiff_t

/* OS definitions and includes */
#define AE_UNKNOWN      0
#define AE_WINDOWS      1
#define AE_POSIX        2
#define AE_LINUX      304

#if AE_OS==AE_WINDOWS
#ifndef _WIN32_WINNT
#define _WIN32_WINNT 0x0501
#endif
#include <windows.h>
#include <process.h>
#elif (AE_OS==AE_POSIX) || (AE_OS==AE_LINUX)
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <sched.h>
#else
#error ALGLIB OS define (AE_OS) is not specified!
#endif

/*************************************************************************
Standard function wrappers for better GLIBC portability
*************************************************************************/
#if defined(_ALGLIB2MKL_GLIBC_WRAPPERS)
#include <math.h>

__asm__(".symver exp,exp@GLIBC_2.2.5");
__asm__(".symver log,log@GLIBC_2.2.5");
__asm__(".symver pow,pow@GLIBC_2.2.5");

double __wrap_exp(double x)
{
    return exp(x);
}

double __wrap_log(double x)
{
    return log(x);
}

double __wrap_pow(double x, double y)
{
    return pow(x, y);
}
#endif

/*
 * This is a watchdog thread which periodically calls mkl_free_buffers()
 * to make sure that no unnecessary memory consumption growth is observed.
 */
static ae_int_t alglib2mkl_initialized = 0;
static ae_int_t watchdog_interval  = 1000;
static ae_int_t watchdog_mem_limit = 25000000;
#if AE_OS==AE_WINDOWS
void watchdog_loop(void *p)
#elif (AE_OS==AE_POSIX) || (AE_OS==AE_LINUX)
void* watchdog_loop(void *p)
#else
#error ALGLIB OS define (AE_OS) is not specified!
#endif
{
    int dummy;
    for(;;)
    {
        /*
         * Sleep for a few seconds
         */
#if AE_OS==AE_WINDOWS
        Sleep((DWORD)watchdog_interval);
#elif (AE_OS==AE_POSIX) || (AE_OS==AE_LINUX)
        {
            struct timespec sleeptime;
            sleeptime.tv_sec  = watchdog_interval/1000;
            sleeptime.tv_nsec = (watchdog_interval%1000)*1000000;
            nanosleep(&sleeptime,NULL);
        }
#else
#error ALGLIB OS define (AE_OS) is not specified!
#endif

        /*
         * Check memory usage; skip buffer deallocation if memory consumption is low.
         */
         if( mkl_mem_stat(&dummy)<watchdog_mem_limit )
             continue;
     
        /*
         * Free buffers
         */
        mkl_free_buffers();
    }
}


/*
 * This function is used under POSIX/Linux to call pthread_create() only once
 */
#if (AE_OS==AE_POSIX) || (AE_OS==AE_LINUX)
static pthread_t watchdog_thread;
void posix_watchdog_primer()
{
    pthread_create(&watchdog_thread, NULL, watchdog_loop, NULL);
    alglib2mkl_initialized = 1;
}
#endif
 

/*
 * This function is used to initialize ALGLIB-MKL interface
 * (in particular, it starts watchdog thread)
 */
#if (AE_OS==AE_POSIX) || (AE_OS==AE_LINUX)
static pthread_once_t watchdog_once_control = PTHREAD_ONCE_INIT;
#endif
void _alglib2mkl_init()
{
    /*
     * Skip if already initialized (OS-agnostic code)
     */
    if( alglib2mkl_initialized )
        return;
    
    /*
     * OS-specific code: try to initialize.
     *
     * Current version of ALGLIB-MKL wrapper allows us to call MKL functions
     * before initialization finalizes, so the only thing we have to do is
     * to make sure that watchdog thread is started only once.
     */
#if AE_OS==AE_WINDOWS
    if( InterlockedCompareExchangePointer((void*)(&alglib2mkl_initialized), (void*)1, (void*)0)!=(void*)0 )
        return;
    _beginthread(watchdog_loop, 0, NULL);
#elif (AE_OS==AE_POSIX) || (AE_OS==AE_LINUX)
    pthread_once(&watchdog_once_control, posix_watchdog_primer);
#else
#error ALGLIB OS define (AE_OS) is not specified!
#endif
}

/*
 * This function is used to perform automatic inplace conversion
 * from ae_int_to to lapack_int:
 * * under LP64 function performs inplace conversion from 64-bit ALGLIB int to 32-bit LAPACK int
 * * under 32-bit environments or ILP64 nothing is done (no conversion is required)
 */
void alglibint_to_lapackint(void *ptr, ptrdiff_t n, char abort_on_failure)
{
    /* volatile is used to trick compiler, to prevent warnings about 'expression is always true' */
    volatile int sz_alglib = 0;
    volatile int sz_lapack = 0;
    
    sz_alglib = sizeof(ptrdiff_t);
    sz_lapack = sizeof(lapack_int);
    
    /* no conversion is required */
    if( sz_alglib==sz_lapack )
        return;
    
    /* ALGLIB int is wider */
    if( sz_alglib==8 && sz_lapack==4 )
    {
        ptrdiff_t i;
        volatile ptrdiff_t  *p_alglib = (ptrdiff_t*)ptr;
        volatile lapack_int *p_lapack = (lapack_int*)ptr;
        
        for(i=0; i<n; i++)
        {
            volatile lapack_int v;
            v = p_alglib[i];
            if( abort_on_failure && v!=p_alglib[i] )
                abort();
            p_lapack[i] = v;
        }
        
        return;
    }
    
    /* unexpected code branch */
    abort();
}


/*
 * This function is used to perform automatic inplace conversion
 * from lapack_int to ae_int_to:
 * * under LP64 function performs inplace conversion from 32-bit LAPACK int to 64-bit ALGLIB int
 * * under 32-bit environments or ILP64 nothing is done (no conversion is required)
 */
void lapackint_to_alglibint(void *ptr, ptrdiff_t n)
{
    /* volatile is used to trick compiler, to prevent warnings about 'expression is always true' */
    volatile int sz_alglib = 0;
    volatile int sz_lapack = 0;
    
    sz_alglib = sizeof(ptrdiff_t);
    sz_lapack = sizeof(lapack_int);
    
    /* no conversion is required */
    if( sz_alglib==sz_lapack )
        return;
    
    /* ALGLIB int is wider */
    if( sz_alglib==8 && sz_lapack==4 )
    {
        ptrdiff_t i;
        volatile ptrdiff_t  *p_alglib = (ptrdiff_t*)ptr;
        volatile lapack_int *p_lapack = (lapack_int*)ptr;
        
        for(i=n-1; i>=0; i--)
        {
            volatile ptrdiff_t v;
            v = p_lapack[i];
            p_alglib[i] = v;
        }
        
        return;
    }
    
    /* unexpected code branch */
    abort();
}

ae_bool _alglib2mkl_disable_fast_mm()
{
    _alglib2mkl_init();
    return mkl_disable_fast_mm();
}

MKL_INT64 _alglib2mkl_memstat()
{
    int dummy;
    _alglib2mkl_init();
    return mkl_mem_stat(&dummy);
}

ae_bool _alglib2mkl_rmatrixtrsv(
     ae_int_t n,
     const double *a_ptr,
     ae_int_t a_stride,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t opa,
     double *x_ptr,
     ae_int_t x_stride)
{
    int prev_nthreads;
    CBLAS_UPLO      uplo  = isupper? CblasUpper : CblasLower;
    CBLAS_DIAG      diag  = isunit ? CblasUnit : CblasNonUnit;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    
    _alglib2mkl_init();
    if( opa==1 )
        trans = CblasTrans;
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dtrsv(CblasRowMajor, uplo, trans, diag, n, a_ptr, a_stride, x_ptr, x_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_cmatrixgerc(
     ptrdiff_t m,
     ptrdiff_t n,
     void* a_ptr,
     ptrdiff_t a_stride,
     const void *alpha,
     const void *u_ptr,
     ptrdiff_t u_stride,
     const void *v_ptr,
     ptrdiff_t v_stride)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_zgerc(CblasRowMajor, m, n, alpha, u_ptr, u_stride, v_ptr, v_stride, a_ptr, a_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_cmatrixgeru(
     ptrdiff_t m,
     ptrdiff_t n,
     void* a_ptr,
     ptrdiff_t a_stride,
     const void *alpha,
     const void *u_ptr,
     ptrdiff_t u_stride,
     const void *v_ptr,
     ptrdiff_t v_stride)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_zgeru(CblasRowMajor, m, n, alpha, u_ptr, u_stride, v_ptr, v_stride, a_ptr, a_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixger(
     ptrdiff_t m,
     ptrdiff_t n,
     double* a_ptr,
     ptrdiff_t a_stride,
     double alpha,
     const double *u_ptr,
     ptrdiff_t u_stride,
     const double *v_ptr,
     ptrdiff_t v_stride)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dger(CblasRowMajor, m, n, alpha, u_ptr, u_stride, v_ptr, v_stride, a_ptr, a_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_cmatrixgemv(
     ae_int_t m,
     ae_int_t n,
     const void *alpha,
     const void *a_ptr,
     ae_int_t a_stride,
     ae_int_t optypea,
     const void *x_ptr,
     ae_int_t x_stride,
     const void *beta,
     void *y_ptr,
     ae_int_t y_stride)
{
    int prev_nthreads;
    CBLAS_TRANSPOSE trans = CblasNoTrans;
    
    _alglib2mkl_init();
    if( optypea==1 )
        trans = CblasTrans;
    if( optypea==2 )
        trans = CblasConjTrans;
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_zgemv(CblasRowMajor, trans, m, n, alpha, a_ptr, a_stride, x_ptr, x_stride, beta, y_ptr, y_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixgemv(
     ae_int_t m,
     ae_int_t n,
     double alpha,
     const double *a_ptr,
     ae_int_t a_stride,
     ae_int_t optypea,
     const double *x_ptr,
     ae_int_t x_stride,
     double beta,
     double *y_ptr,
     ae_int_t y_stride)
{
    int prev_nthreads;
    CBLAS_TRANSPOSE trans = optypea==0 ? CblasNoTrans : CblasTrans;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dgemv(CblasRowMajor, trans, m, n, alpha, a_ptr, a_stride, x_ptr, x_stride, beta, y_ptr, y_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixsymv(
     ae_int_t n,
     double alpha,
     const double *a_ptr,
     ae_int_t a_stride,
     ae_bool isupper,
     const double *x_ptr,
     ae_int_t x_stride,
     double beta,
     double *y_ptr,
     ae_int_t y_stride)
{
    int prev_nthreads;
    CBLAS_UPLO uplo = isupper? CblasUpper : CblasLower;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dsymv(CblasRowMajor, uplo, n, alpha, a_ptr, a_stride, x_ptr, x_stride, beta, y_ptr, y_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixsyrkmkl(
     ae_int_t n,
     ae_int_t k,
     double alpha,
     double *a_ptr,
     ae_int_t a_stride,
     ae_int_t optypea,
     double beta,
     double *c_ptr,
     ae_int_t c_stride,
     ae_bool isupper)
{
    int prev_nthreads;
    CBLAS_UPLO uplo = isupper ? CblasUpper : CblasLower;
    CBLAS_TRANSPOSE trans = optypea==0 ? CblasNoTrans : CblasTrans;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dsyrk(
        CblasRowMajor, uplo, trans, n, k,
        alpha, a_ptr, a_stride,
        beta,  c_ptr, c_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_cmatrixherkmkl(
     ae_int_t n,
     ae_int_t k,
     double alpha,
     void *a_ptr,
     ae_int_t a_stride,
     ae_int_t optypea,
     double beta,
     void *c_ptr,
     ae_int_t c_stride,
     ae_bool isupper)
{
    int prev_nthreads;
    CBLAS_UPLO uplo = isupper ? CblasUpper : CblasLower;
    CBLAS_TRANSPOSE trans = optypea==0 ? CblasNoTrans : CblasConjTrans;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_zherk(
        CblasRowMajor, uplo, trans, n, k,
        alpha, a_ptr, a_stride,
        beta,  c_ptr, c_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixgemmmkl(
     ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     double *a_ptr,
     ae_int_t a_stride,
     ae_int_t optypea,
     double *b_ptr,
     ae_int_t b_stride,
     ae_int_t optypeb,
     double beta,
     double *c_ptr,
     ae_int_t c_stride)
{
    CBLAS_TRANSPOSE trans_a = optypea==0 ? CblasNoTrans : CblasTrans;
    CBLAS_TRANSPOSE trans_b = optypeb==0 ? CblasNoTrans : CblasTrans;
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dgemm(
        CblasRowMajor, trans_a, trans_b,
        m, n, k,
        alpha,
        a_ptr, a_stride,
        b_ptr, b_stride,
        beta,
        c_ptr, c_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_cmatrixgemmmkl(
     ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     void* alpha,
     void* a_ptr,
     ae_int_t a_stride,
     ae_int_t optypea,
     void* b_ptr,
     ae_int_t b_stride,
     ae_int_t optypeb,
     void* beta,
     void* c_ptr,
     ae_int_t c_stride)
{
    CBLAS_TRANSPOSE trans_a = optypea==0 ? CblasNoTrans : (optypea==1 ? CblasTrans : CblasConjTrans);
    CBLAS_TRANSPOSE trans_b = optypeb==0 ? CblasNoTrans : (optypeb==1 ? CblasTrans : CblasConjTrans);
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_zgemm(
        CblasRowMajor, trans_a, trans_b,
        m, n, k,
        alpha,
        a_ptr, a_stride,
        b_ptr, b_stride,
        beta,
        c_ptr, c_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixtrsmmkl(
     char is_rightside,
     ptrdiff_t m,
     ptrdiff_t n,
     double* a,
     ptrdiff_t a_stride,
     char isupper,
     char isunit,
     ptrdiff_t optype,
     double* x,
     ptrdiff_t x_stride)
{
    CBLAS_SIDE cb_side = is_rightside ? CblasRight : CblasLeft;
    CBLAS_UPLO cb_uplo = isupper ? CblasUpper : CblasLower;
    CBLAS_DIAG cb_diag = isunit ? CblasUnit : CblasNonUnit;
    CBLAS_TRANSPOSE cb_trans_a = optype==0 ? CblasNoTrans : CblasTrans;
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_dtrsm(
        CblasRowMajor,
        cb_side, cb_uplo, cb_trans_a, cb_diag,
        m, n,
        1.0,
        a, a_stride,
        x, x_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_cmatrixtrsmmkl(
     char is_rightside,
     ptrdiff_t m,
     ptrdiff_t n,
     void* a,
     ptrdiff_t a_stride,
     char isupper,
     char isunit,
     ptrdiff_t optype,
     void* x,
     ptrdiff_t x_stride)
{
    CBLAS_SIDE cb_side = is_rightside ? CblasRight : CblasLeft;
    CBLAS_UPLO cb_uplo = isupper ? CblasUpper : CblasLower;
    CBLAS_DIAG cb_diag = isunit ? CblasUnit : CblasNonUnit;
    CBLAS_TRANSPOSE cb_trans_a = optype==0 ? CblasNoTrans : (optype==1 ? CblasTrans : CblasConjTrans);
    int prev_nthreads;
    double alpha[2] = { 1, 0 };
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    cblas_ztrsm(
        CblasRowMajor,
        cb_side, cb_uplo, cb_trans_a, cb_diag,
        m, n,
        (void*)alpha,
        a, a_stride,
        x, x_stride);
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_spdmatrixcholeskymkl(
     double* a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     ae_bool  isupper,
     ae_bool* cholresult)
{
    char uplo = isupper ? 'U' : 'L';
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    *cholresult = LAPACKE_dpotrf(LAPACK_ROW_MAJOR, uplo, n, a, a_stride)==0;
    mkl_set_num_threads_local(prev_nthreads);
    return ae_true;
}

ae_bool _alglib2mkl_rmatrixplumkl(
     double* a,
     ptrdiff_t a_stride,
     ptrdiff_t m,
     ptrdiff_t n,
     ptrdiff_t offs,
     ptrdiff_t *pivots)
{
    ptrdiff_t minmn, i;
    int prev_nthreads;
    
    _alglib2mkl_init();
    
    /* call LAPACK code */
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, a, a_stride, (lapack_int*)pivots);
    mkl_set_num_threads_local(prev_nthreads);
    
    /* convert from LAPACK numeration back to ALGLIB */
    minmn = m<n ? m : n;
    lapackint_to_alglibint(pivots, minmn);
    for(i=0; i<minmn; i++)
        pivots[i] = (pivots[i]-1)+offs;
    return ae_true;
}

char _alglib2mkl_rmatrixbdmkl(
     double    *a,
     ptrdiff_t a_stride,
     ptrdiff_t m,
     ptrdiff_t n,
     double    *d,
     double    *e,
     double    *tauq,
     double    *taup)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dgebrd(LAPACK_ROW_MAJOR, m, n, a, a_stride, d, e, tauq, taup);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_rmatrixbdmultiplybymkl(
     double *qp,
     ptrdiff_t qp_stride,
     ptrdiff_t m,
     ptrdiff_t n,
     double *tauq,
     double *taup,
     double *z,
     ptrdiff_t z_stride,
     ptrdiff_t zrows,
     ptrdiff_t zcolumns,
     char byq,
     char fromtheright,
     char dotranspose)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dormbr(
        LAPACK_ROW_MAJOR,
        byq ? 'Q' : 'P',
        fromtheright ? 'R' : 'L',
        dotranspose  ? 'T' : 'N',
        zrows, zcolumns,
        byq ? n : m,
        qp, qp_stride,
        byq ? tauq : taup,
        z, z_stride);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_rmatrixhessenbergmkl(
     double *a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     double* tau)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dgehrd(LAPACK_ROW_MAJOR, n, 1, n, a, a_stride, tau);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_rmatrixhessenbergunpackqmkl(
     double    *a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     double   *tau,
     double    *q,
     ptrdiff_t q_stride)
{
    int prev_nthreads, i, j;
    
    _alglib2mkl_init();
    for(i=0; i<n; i++)
    {
        double *p_arow = a+a_stride*i;
        double *p_qrow = q+q_stride*i;
        for(j=0; j<n; j++)
            p_qrow[j] = p_arow[j];
    }
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dorghr(LAPACK_ROW_MAJOR, n, 1, n, q, q_stride, tau);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_smatrixtdmkl(
     double *a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     char isupper,
     double* tau,
     double* d,
     double* e)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dsytrd(
        LAPACK_ROW_MAJOR,
        isupper ? 'U' : 'L',
        n,
        a, a_stride,
        d, e, tau);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_smatrixtdunpackqmkl(
     double* a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     char isupper,
     double* tau,
     double* q,
     ptrdiff_t q_stride)
{
    int prev_nthreads, i, j;
    
    _alglib2mkl_init();
    for(i=0; i<n; i++)
    {
        double *p_arow = a+a_stride*i;
        double *p_qrow = q+q_stride*i;
        for(j=0; j<n; j++)
            p_qrow[j] = p_arow[j];
    }
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_dorgtr(
        LAPACK_ROW_MAJOR,
        isupper ? 'U' : 'L',
        n,
        q, q_stride,
        tau);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_hmatrixtdmkl(
     void *a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     char isupper,
     void* tau,
     double* d,
     double* e)
{
    int prev_nthreads;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_zhetrd(
        LAPACK_ROW_MAJOR,
        isupper ? 'U' : 'L',
        n,
        (lapack_complex_double*)a, a_stride,
        d, e, (lapack_complex_double*)tau);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_hmatrixtdunpackqmkl(
     void* a,
     ptrdiff_t a_stride,
     ptrdiff_t n,
     char isupper,
     void* tau,
     void* q,
     ptrdiff_t q_stride)
{
    int prev_nthreads, i, j;
    
    _alglib2mkl_init();
    for(i=0; i<n; i++)
    {
        double *p_arow = (double*)a+2*a_stride*i;
        double *p_qrow = (double*)q+2*q_stride*i;
        for(j=0; j<2*n; j++)
            p_qrow[j] = p_arow[j];
    }
    prev_nthreads = mkl_set_num_threads_local(1);
    LAPACKE_zungtr(
        LAPACK_ROW_MAJOR,
        isupper ? 'U' : 'L',
        n,
        (lapack_complex_double*)q, q_stride,
        (lapack_complex_double*)tau);
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_rmatrixbdsvdmkl(
     double* d,
     double* e,
     ptrdiff_t n,
     char isupper,
     double* u,
     ptrdiff_t u_stride,
     ptrdiff_t nru,
     double* c,
     ptrdiff_t c_stride,
     ptrdiff_t ncc,
     double* vt,
     ptrdiff_t vt_stride,
     ptrdiff_t ncvt,
     char* svdresult)
{
    int prev_nthreads;
    lapack_int result;
    
    _alglib2mkl_init();
    prev_nthreads = mkl_set_num_threads_local(1);
    result = LAPACKE_dbdsqr(
        LAPACK_ROW_MAJOR,
        isupper ? 'U' : 'L',
        n, (lapack_int)ncvt, (lapack_int)nru, (lapack_int)ncc,
        d, e, vt, (lapack_int)vt_stride, u, (lapack_int)u_stride, c, (lapack_int)c_stride);
    *svdresult = result==0;
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_rmatrixinternalschurdecompositionmkl(
     double* h,
     ptrdiff_t h_stride,
     ptrdiff_t n,
     ptrdiff_t tneeded,
     ptrdiff_t zneeded,
     double* wr,
     double* wi,
     double* z,
     ptrdiff_t z_stride,
     ptrdiff_t* info)
{
    int prev_nthreads;
    lapack_int result;
    char compz;
    
    _alglib2mkl_init();
    compz = 'N';
    if( zneeded==1 )
        compz = 'V';
    if( zneeded==2 )
        compz = 'I';
    prev_nthreads = mkl_set_num_threads_local(1);
    result = LAPACKE_dhseqr(
        LAPACK_ROW_MAJOR,
        tneeded ? 'S' : 'E',
        compz,
        (lapack_int)n, 1, (lapack_int)n,
        h, (lapack_int)h_stride,
        wr, wi,
        z, (lapack_int)z_stride);
    *info = result;
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_rmatrixinternaltrevcmkl(
     double* t,
     ptrdiff_t t_stride,
     ptrdiff_t n,
     ptrdiff_t side,
     ptrdiff_t howmny,
     double* vl,
     ptrdiff_t vl_stride,
     double* vr,
     ptrdiff_t vr_stride,
     ptrdiff_t* m,
     ptrdiff_t* info)
{
    int prev_nthreads;
    lapack_int loc_m;
    char c_side, c_howmny;

    _alglib2mkl_init();
    if( howmny==3 ) /* not supported! */
        return ae_false;
    c_side = 'B';
    if( side==1 )
        c_side = 'R';
    if( side==2 )
        c_side = 'L';
    c_howmny = howmny==1 ? 'B' : 'A';
    prev_nthreads = mkl_set_num_threads_local(1);
    *info = LAPACKE_dtrevc(
        LAPACK_ROW_MAJOR,
        c_side,
        c_howmny,
        NULL,
        n, 
        t, t_stride,
        vl, vl_stride,
        vr, vr_stride,
        n,
        &loc_m);
    *m = loc_m;
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_smatrixtdevdmkl(
     double* d,
     double* e,
     ptrdiff_t n,
     ptrdiff_t zneeded,
     double* z,
     ptrdiff_t z_stride,
     char* evdresult)
{
    lapack_int r;
    int prev_nthreads;
    char compz;
    
    _alglib2mkl_init();
    
    /* only zneeded=0,1,2 are supported */
    if( zneeded==3 )
        return ae_false;
        
    /* general case */
    prev_nthreads = mkl_set_num_threads_local(1);
    compz = 'N';
    if( zneeded==1 )
        compz = 'V';
    if( zneeded==2 )
        compz = 'I';
    r = LAPACKE_dsteqr(
        LAPACK_ROW_MAJOR,
        compz,
        (lapack_int)n, d, e,
        z, (lapack_int)z_stride);
    *evdresult = r==0;
    mkl_set_num_threads_local(prev_nthreads);
    
    return ae_true;
}

char _alglib2mkl_sparsegemvcrsmkl(
     ptrdiff_t opa,
     ptrdiff_t arows,
     ptrdiff_t acols,
     double alpha,
     const double* vals,
     const ptrdiff_t* cidx,
     const ptrdiff_t* ridx,
     const double* x,
     double beta,
     double* y)
{
    char matdescra[4];
    
    /*
     * LP64 is not supported
     */
    if( sizeof(MKL_INT)!=sizeof(ptrdiff_t) )
        return ae_false;

    /*
     * Perform initialization
     */
    _alglib2mkl_init();
    
    /*
     * Call MKL
     */
    matdescra[0] = 'G';
    matdescra[1] = '?';
    matdescra[2] = '?';
    matdescra[3] = 'C';
    mkl_dcsrmv(opa==0 ? "N" : "T", (MKL_INT*)(&arows), (MKL_INT*)(&acols),
        &alpha,
        matdescra, vals, (const MKL_INT*)cidx , (const MKL_INT *)ridx, (const MKL_INT *)(ridx+1),
        x,
        &beta,
        y);
    return ae_true;
}
