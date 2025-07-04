/*************************************************************************
ALGLIB 3.19.0 (source code generated 2022-06-05)
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This  program  is  a  Commercial Edition  of the ALGLIB package  licensed
to Salman (Licensee), agreement ID is AGR-20211124-1

As of 2022-10-12:
* DEV-1 license plan is purchased (1 developer)
* following developers are registered for this license:
  * umut.salman@lspm.cnrs.fr

========================== GENERAL INFORMATION ==========================

1. Only Licensee and  its Sublicensees  can  use/distribute  it according
to ALGLIB License Agreement (see below) between Licensor and Licensee.

2. All developers working  for  Licensee  should  register  themselves at
alglib.net.

3. This  source  code  may  contain  modifications  made  by  Licensee or
Sublicensees  which  fall under the terms of the ALGLIB License Agreement
too.

4. Text  below  is  an excerpt from complete ALGLIB License Agreement,  a
part which governs usage and redistribution rights granted  to  Licensee.
See agreement-v8uk.pdf at the root of  ALGLIB distribution for a complete
text of the license agreement.

================ ALGLIB LICENSE AGREEMENT ( APPENDIX A ) ================

DEFINITIONS:
* "ALGLIB" – software delivered by Licensor  to  Licensee  under  present
  Agreement. ALGLIB may include  Binary  Components  (delivered  only  in
  binary form) and Source  Code  Components  (with  optional  precompiled
  binary form).  ALGLIB  includes  integrated  third  party  software  as
  defined by clause 1 of the Appendix B which is considered as a part  of
  ALGLIB for the purposes of the present Agreement. Optional third  party
  software (as defined by clause 2 of Appendix B) is not considered as  a
  part of ALGLIB (even when such software is bundled with ALGLIB) and has
  its own licensing terms.
* "Application"  -  program  developed  by  Licensee  (either  standalone
  application or software development library) which includes  ALGLIB  as
  one of its parts .
* "Sublicensee"  -  any  party  (including  resellers)   which   receives
  Application from Licensee or another Sublicensee.
* "Application License Agreement"  -  agreement  which   governs   usage/
  redistribution of the Application.
  
LICENSE GRANT:
Subject to the License Restrictions below, Licensor  grants  to  Licensee
the following non-exclusive royalty-free licenses:
A. To modify Source Code Components of ALGLIB and to use modified version
   on the terms of this Agreement.
B. To  develop  Applications  which  use  ALGLIB  and  to distribute such
   Applications in Binary and/or Source Code forms,  with  ALGLIB  either
   statically or dynamically linked. This right is granted provided that:
   * distribution of Source Code forms of Application/ALGLIB is performed
     subject to additional conditions set by clause H (this clause is not
     applied to binary-only distribution)
   * such Applications add significant  primary  functionality  different
     from that of the ALGLIB.
   * such Applications do not expose ALGLIB API (application  programming
     interface) either directly or indirectly
   * Sublicensee  has  no   right   to  use  ALGLIB  except  as  part  of
     the Application
   * any  subsequent  redistribution   respects    conditions    of   the
     present Agreement
   * all Licensee’s developers using ALGLIB should register at  company's
     account at www.alglib.net
C. To use Resellers for distribution of the  Application  (in  Binary  or
   Source Code forms), provided that the only activity Reseller  performs
   with Application is redistribution.
   
LICENSE RESTRICTIONS:
D. Licensee/Sublicensee may NOT use, copy or distribute ALGLIB except  as
   provided in this Agreement.
D2. Licensee/Sublicensee may NOT rent or lease ALGLIB to any third party.
E. Licensee/Sublicensee may NOT disassemble, reverse engineer, decompile,
   modify Binary Components of ALGLIB or compiled forms  of  Source  Code
   components.
F. Licensee/Sublicensee  may  NOT  remove  any  copyright notice from the
   Source Code / Binary Components.
G. Licensee/Sublicensee may NOT  disable/remove  code  which  checks  for
   presence of license keys (if such code is included in ALGLIB) from the
   Source Code / Binary Components.
H. Distribution of  Source  Code  forms  of  Application/ALGLIB  must  be
   performed subject to additional conditions:
   * Source Code Components of ALGLIB are distributed only as part of the
     Application. They are not  publicly  distributed.  Sublicensee  must
     explicitly accept Application License Agreement in order  to  access
     ALGLIB source code.
   * Sublicensee has no right to redistribute Application/ALGLIB (in  any
     form, Binary or Source Code), unless Sublicensee is Reseller who  is
     fully compliant with conditions set by clause C.
   * Sublicensee has no right to modify ALGLIB Source  Code,  except  for
     the purpose of fixing bugs
   * Sublicensee has no right to workaround "use ALGLIB only as  part  of
     the Application" limitation by sequentially modifying Application in
     a way which effectively creates new program with different  purpose.
     Application   License  Agreement  may  (a)  explicitly  forbid  such
     modifications, or (b) allow only limited set of "safe" modifications
     (developing plugins, fixing bugs, modifying only specific  parts  of
     the Application).
     
COPYRIGHT:
Title to the ALGLIB and all copies  thereof  remain with Licensor. ALGLIB
is copyrighted and is protected by international copyright laws. You will
not remove any copyright notice from  the  ALGLIB  files.  You  agree  to
prevent any unauthorized copying of ALGLIB. Except as expressly  provided
herein, Licensor does not grant any express or implied right to you under
Licensor patents, copyrights, trademarks, or trade secret information.
>>> END OF LICENSE >>>
*************************************************************************/
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "stdafx.h"
#include "diffequations.h"

// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC) && !defined(AE_ALL_WARNINGS)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4611)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{

#if defined(AE_COMPILE_ODESOLVER) || !defined(AE_PARTIAL_BUILD)

#endif

#if defined(AE_COMPILE_ODESOLVER) || !defined(AE_PARTIAL_BUILD)
/*************************************************************************

*************************************************************************/
_odesolverstate_owner::_odesolverstate_owner()
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _state;
    
    alglib_impl::ae_state_init(&_state);
    if( setjmp(_break_jump) )
    {
        if( p_struct!=NULL )
        {
            alglib_impl::_odesolverstate_destroy(p_struct);
            alglib_impl::ae_free(p_struct);
        }
        p_struct = NULL;
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_state.error_msg);
        return;
#endif
    }
    alglib_impl::ae_state_set_break_jump(&_state, &_break_jump);
    p_struct = NULL;
    p_struct = (alglib_impl::odesolverstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::odesolverstate), &_state);
    memset(p_struct, 0, sizeof(alglib_impl::odesolverstate));
    alglib_impl::_odesolverstate_init(p_struct, &_state, ae_false);
    ae_state_clear(&_state);
}

_odesolverstate_owner::_odesolverstate_owner(const _odesolverstate_owner &rhs)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _state;
    
    alglib_impl::ae_state_init(&_state);
    if( setjmp(_break_jump) )
    {
        if( p_struct!=NULL )
        {
            alglib_impl::_odesolverstate_destroy(p_struct);
            alglib_impl::ae_free(p_struct);
        }
        p_struct = NULL;
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_state.error_msg);
        return;
#endif
    }
    alglib_impl::ae_state_set_break_jump(&_state, &_break_jump);
    p_struct = NULL;
    alglib_impl::ae_assert(rhs.p_struct!=NULL, "ALGLIB: odesolverstate copy constructor failure (source is not initialized)", &_state);
    p_struct = (alglib_impl::odesolverstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::odesolverstate), &_state);
    memset(p_struct, 0, sizeof(alglib_impl::odesolverstate));
    alglib_impl::_odesolverstate_init_copy(p_struct, const_cast<alglib_impl::odesolverstate*>(rhs.p_struct), &_state, ae_false);
    ae_state_clear(&_state);
}

_odesolverstate_owner& _odesolverstate_owner::operator=(const _odesolverstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    jmp_buf _break_jump;
    alglib_impl::ae_state _state;
    
    alglib_impl::ae_state_init(&_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_state.error_msg);
        return *this;
#endif
    }
    alglib_impl::ae_state_set_break_jump(&_state, &_break_jump);
    alglib_impl::ae_assert(p_struct!=NULL, "ALGLIB: odesolverstate assignment constructor failure (destination is not initialized)", &_state);
    alglib_impl::ae_assert(rhs.p_struct!=NULL, "ALGLIB: odesolverstate assignment constructor failure (source is not initialized)", &_state);
    alglib_impl::_odesolverstate_destroy(p_struct);
    memset(p_struct, 0, sizeof(alglib_impl::odesolverstate));
    alglib_impl::_odesolverstate_init_copy(p_struct, const_cast<alglib_impl::odesolverstate*>(rhs.p_struct), &_state, ae_false);
    ae_state_clear(&_state);
    return *this;
}

_odesolverstate_owner::~_odesolverstate_owner()
{
    if( p_struct!=NULL )
    {
        alglib_impl::_odesolverstate_destroy(p_struct);
        ae_free(p_struct);
    }
}

alglib_impl::odesolverstate* _odesolverstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::odesolverstate* _odesolverstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::odesolverstate*>(p_struct);
}
odesolverstate::odesolverstate() : _odesolverstate_owner() ,needdy(p_struct->needdy),y(&p_struct->y),dy(&p_struct->dy),x(p_struct->x)
{
}

odesolverstate::odesolverstate(const odesolverstate &rhs):_odesolverstate_owner(rhs) ,needdy(p_struct->needdy),y(&p_struct->y),dy(&p_struct->dy),x(p_struct->x)
{
}

odesolverstate& odesolverstate::operator=(const odesolverstate &rhs)
{
    if( this==&rhs )
        return *this;
    _odesolverstate_owner::operator=(rhs);
    return *this;
}

odesolverstate::~odesolverstate()
{
}


/*************************************************************************

*************************************************************************/
_odesolverreport_owner::_odesolverreport_owner()
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _state;
    
    alglib_impl::ae_state_init(&_state);
    if( setjmp(_break_jump) )
    {
        if( p_struct!=NULL )
        {
            alglib_impl::_odesolverreport_destroy(p_struct);
            alglib_impl::ae_free(p_struct);
        }
        p_struct = NULL;
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_state.error_msg);
        return;
#endif
    }
    alglib_impl::ae_state_set_break_jump(&_state, &_break_jump);
    p_struct = NULL;
    p_struct = (alglib_impl::odesolverreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::odesolverreport), &_state);
    memset(p_struct, 0, sizeof(alglib_impl::odesolverreport));
    alglib_impl::_odesolverreport_init(p_struct, &_state, ae_false);
    ae_state_clear(&_state);
}

_odesolverreport_owner::_odesolverreport_owner(const _odesolverreport_owner &rhs)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _state;
    
    alglib_impl::ae_state_init(&_state);
    if( setjmp(_break_jump) )
    {
        if( p_struct!=NULL )
        {
            alglib_impl::_odesolverreport_destroy(p_struct);
            alglib_impl::ae_free(p_struct);
        }
        p_struct = NULL;
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_state.error_msg);
        return;
#endif
    }
    alglib_impl::ae_state_set_break_jump(&_state, &_break_jump);
    p_struct = NULL;
    alglib_impl::ae_assert(rhs.p_struct!=NULL, "ALGLIB: odesolverreport copy constructor failure (source is not initialized)", &_state);
    p_struct = (alglib_impl::odesolverreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::odesolverreport), &_state);
    memset(p_struct, 0, sizeof(alglib_impl::odesolverreport));
    alglib_impl::_odesolverreport_init_copy(p_struct, const_cast<alglib_impl::odesolverreport*>(rhs.p_struct), &_state, ae_false);
    ae_state_clear(&_state);
}

_odesolverreport_owner& _odesolverreport_owner::operator=(const _odesolverreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    jmp_buf _break_jump;
    alglib_impl::ae_state _state;
    
    alglib_impl::ae_state_init(&_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_state.error_msg);
        return *this;
#endif
    }
    alglib_impl::ae_state_set_break_jump(&_state, &_break_jump);
    alglib_impl::ae_assert(p_struct!=NULL, "ALGLIB: odesolverreport assignment constructor failure (destination is not initialized)", &_state);
    alglib_impl::ae_assert(rhs.p_struct!=NULL, "ALGLIB: odesolverreport assignment constructor failure (source is not initialized)", &_state);
    alglib_impl::_odesolverreport_destroy(p_struct);
    memset(p_struct, 0, sizeof(alglib_impl::odesolverreport));
    alglib_impl::_odesolverreport_init_copy(p_struct, const_cast<alglib_impl::odesolverreport*>(rhs.p_struct), &_state, ae_false);
    ae_state_clear(&_state);
    return *this;
}

_odesolverreport_owner::~_odesolverreport_owner()
{
    if( p_struct!=NULL )
    {
        alglib_impl::_odesolverreport_destroy(p_struct);
        ae_free(p_struct);
    }
}

alglib_impl::odesolverreport* _odesolverreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::odesolverreport* _odesolverreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::odesolverreport*>(p_struct);
}
odesolverreport::odesolverreport() : _odesolverreport_owner() ,nfev(p_struct->nfev),terminationtype(p_struct->terminationtype)
{
}

odesolverreport::odesolverreport(const odesolverreport &rhs):_odesolverreport_owner(rhs) ,nfev(p_struct->nfev),terminationtype(p_struct->terminationtype)
{
}

odesolverreport& odesolverreport::operator=(const odesolverreport &rhs)
{
    if( this==&rhs )
        return *this;
    _odesolverreport_owner::operator=(rhs);
    return *this;
}

odesolverreport::~odesolverreport()
{
}

/*************************************************************************
Cash-Karp adaptive ODE solver.

This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
(here Y may be single variable or vector of N variables).

INPUT PARAMETERS:
    Y       -   initial conditions, array[0..N-1].
                contains values of Y[] at X[0]
    N       -   system size
    X       -   points at which Y should be tabulated, array[0..M-1]
                integrations starts at X[0], ends at X[M-1],  intermediate
                values at X[i] are returned too.
                SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!
    M       -   number of intermediate points + first point + last point:
                * M>2 means that you need both Y(X[M-1]) and M-2 values at
                  intermediate points
                * M=2 means that you want just to integrate from  X[0]  to
                  X[1] and don't interested in intermediate values.
                * M=1 means that you don't want to integrate :)
                  it is degenerate case, but it will be handled correctly.
                * M<1 means error
    Eps     -   tolerance (absolute/relative error on each  step  will  be
                less than Eps). When passing:
                * Eps>0, it means desired ABSOLUTE error
                * Eps<0, it means desired RELATIVE error.  Relative errors
                  are calculated with respect to maximum values of  Y seen
                  so far. Be careful to use this criterion  when  starting
                  from Y[] that are close to zero.
    H       -   initial  step  lenth,  it  will  be adjusted automatically
                after the first  step.  If  H=0,  step  will  be  selected
                automatically  (usualy  it  will  be  equal  to  0.001  of
                min(x[i]-x[j])).

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of OdeSolverIteration. Used for reverse communication.
                This structure should be passed  to the OdeSolverIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
void odesolverrkck(const real_1d_array &y, const ae_int_t n, const real_1d_array &x, const ae_int_t m, const double eps, const double h, odesolverstate &state, const xparams _xparams)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
        return;
#endif
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    if( _xparams.flags!=0x0 )
        ae_state_set_flags(&_alglib_env_state, _xparams.flags);
    alglib_impl::odesolverrkck(const_cast<alglib_impl::ae_vector*>(y.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), m, eps, h, const_cast<alglib_impl::odesolverstate*>(state.c_ptr()), &_alglib_env_state);
    alglib_impl::ae_state_clear(&_alglib_env_state);
    return;
}

/*************************************************************************
Cash-Karp adaptive ODE solver.

This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
(here Y may be single variable or vector of N variables).

INPUT PARAMETERS:
    Y       -   initial conditions, array[0..N-1].
                contains values of Y[] at X[0]
    N       -   system size
    X       -   points at which Y should be tabulated, array[0..M-1]
                integrations starts at X[0], ends at X[M-1],  intermediate
                values at X[i] are returned too.
                SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!
    M       -   number of intermediate points + first point + last point:
                * M>2 means that you need both Y(X[M-1]) and M-2 values at
                  intermediate points
                * M=2 means that you want just to integrate from  X[0]  to
                  X[1] and don't interested in intermediate values.
                * M=1 means that you don't want to integrate :)
                  it is degenerate case, but it will be handled correctly.
                * M<1 means error
    Eps     -   tolerance (absolute/relative error on each  step  will  be
                less than Eps). When passing:
                * Eps>0, it means desired ABSOLUTE error
                * Eps<0, it means desired RELATIVE error.  Relative errors
                  are calculated with respect to maximum values of  Y seen
                  so far. Be careful to use this criterion  when  starting
                  from Y[] that are close to zero.
    H       -   initial  step  lenth,  it  will  be adjusted automatically
                after the first  step.  If  H=0,  step  will  be  selected
                automatically  (usualy  it  will  be  equal  to  0.001  of
                min(x[i]-x[j])).

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of OdeSolverIteration. Used for reverse communication.
                This structure should be passed  to the OdeSolverIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
#if !defined(AE_NO_EXCEPTIONS)
void odesolverrkck(const real_1d_array &y, const real_1d_array &x, const double eps, const double h, odesolverstate &state, const xparams _xparams)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;
    ae_int_t m;

    n = y.length();
    m = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
        _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    if( _xparams.flags!=0x0 )
        ae_state_set_flags(&_alglib_env_state, _xparams.flags);
    alglib_impl::odesolverrkck(const_cast<alglib_impl::ae_vector*>(y.c_ptr()), n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), m, eps, h, const_cast<alglib_impl::odesolverstate*>(state.c_ptr()), &_alglib_env_state);

    alglib_impl::ae_state_clear(&_alglib_env_state);
    return;
}
#endif

/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool odesolveriteration(const odesolverstate &state, const xparams _xparams)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
        return 0;
#endif
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    if( _xparams.flags!=0x0 )
        ae_state_set_flags(&_alglib_env_state, _xparams.flags);
    ae_bool result = alglib_impl::odesolveriteration(const_cast<alglib_impl::odesolverstate*>(state.c_ptr()), &_alglib_env_state);
    alglib_impl::ae_state_clear(&_alglib_env_state);
    return *(reinterpret_cast<bool*>(&result));
}


void odesolversolve(odesolverstate &state,
    void (*diff)(const real_1d_array &y, double x, real_1d_array &dy, void *ptr),
    void *ptr, const xparams _xparams){
    jmp_buf _break_jump;
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
        return;
#endif
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    if( _xparams.flags!=0x0 )
        ae_state_set_flags(&_alglib_env_state, _xparams.flags);
    alglib_impl::ae_assert(diff!=NULL, "ALGLIB: error in 'odesolversolve()' (diff is NULL)", &_alglib_env_state);
    while( alglib_impl::odesolveriteration(state.c_ptr(), &_alglib_env_state) )
    {
        _ALGLIB_CALLBACK_EXCEPTION_GUARD_BEGIN
                if( state.needdy )
                {
                    diff(state.y, state.x, state.dy, ptr);
                    continue;
                }
        goto lbl_no_callback;
        _ALGLIB_CALLBACK_EXCEPTION_GUARD_END
    lbl_no_callback:
        alglib_impl::ae_assert(ae_false, "ALGLIB: unexpected error in 'odesolversolve'", &_alglib_env_state);
    }
    alglib_impl::ae_state_clear(&_alglib_env_state);
}



/*************************************************************************
ODE solver results

Called after OdeSolverIteration returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by OdeSolverIteration).

OUTPUT PARAMETERS:
    M       -   number of tabulated values, M>=1
    XTbl    -   array[0..M-1], values of X
    YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]
    Rep     -   solver report:
                * Rep.TerminationType completetion code:
                    * -2    X is not ordered  by  ascending/descending  or
                            there are non-distinct X[],  i.e.  X[i]=X[i+1]
                    * -1    incorrect parameters were specified
                    *  1    task has been solved
                * Rep.NFEV contains number of function calculations

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
void odesolverresults(const odesolverstate &state, ae_int_t &m, real_1d_array &xtbl, real_2d_array &ytbl, odesolverreport &rep, const xparams _xparams)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
        return;
#endif
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    if( _xparams.flags!=0x0 )
        ae_state_set_flags(&_alglib_env_state, _xparams.flags);
    alglib_impl::odesolverresults(const_cast<alglib_impl::odesolverstate*>(state.c_ptr()), &m, const_cast<alglib_impl::ae_vector*>(xtbl.c_ptr()), const_cast<alglib_impl::ae_matrix*>(ytbl.c_ptr()), const_cast<alglib_impl::odesolverreport*>(rep.c_ptr()), &_alglib_env_state);
    alglib_impl::ae_state_clear(&_alglib_env_state);
    return;
}
#endif
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF COMPUTATIONAL CORE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
#if defined(AE_COMPILE_ODESOLVER) || !defined(AE_PARTIAL_BUILD)
static double odesolver_odesolvermaxgrow = 3.0;
static double odesolver_odesolvermaxshrink = 10.0;
static double odesolver_odesolverguaranteeddecay = 0.9;
static void odesolver_odesolverinit(ae_int_t solvertype,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_int_t m,
     double eps,
     double h,
     odesolverstate* state,
     ae_state *_state);


#endif

#if defined(AE_COMPILE_ODESOLVER) || !defined(AE_PARTIAL_BUILD)


/*************************************************************************
Cash-Karp adaptive ODE solver.

This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
(here Y may be single variable or vector of N variables).

INPUT PARAMETERS:
    Y       -   initial conditions, array[0..N-1].
                contains values of Y[] at X[0]
    N       -   system size
    X       -   points at which Y should be tabulated, array[0..M-1]
                integrations starts at X[0], ends at X[M-1],  intermediate
                values at X[i] are returned too.
                SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!
    M       -   number of intermediate points + first point + last point:
                * M>2 means that you need both Y(X[M-1]) and M-2 values at
                  intermediate points
                * M=2 means that you want just to integrate from  X[0]  to
                  X[1] and don't interested in intermediate values.
                * M=1 means that you don't want to integrate :)
                  it is degenerate case, but it will be handled correctly.
                * M<1 means error
    Eps     -   tolerance (absolute/relative error on each  step  will  be
                less than Eps). When passing:
                * Eps>0, it means desired ABSOLUTE error
                * Eps<0, it means desired RELATIVE error.  Relative errors
                  are calculated with respect to maximum values of  Y seen
                  so far. Be careful to use this criterion  when  starting
                  from Y[] that are close to zero.
    H       -   initial  step  lenth,  it  will  be adjusted automatically
                after the first  step.  If  H=0,  step  will  be  selected
                automatically  (usualy  it  will  be  equal  to  0.001  of
                min(x[i]-x[j])).

OUTPUT PARAMETERS
    State   -   structure which stores algorithm state between  subsequent
                calls of OdeSolverIteration. Used for reverse communication.
                This structure should be passed  to the OdeSolverIteration
                subroutine.

SEE ALSO
    AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.


  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
void odesolverrkck(/* Real    */ ae_vector* y,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_int_t m,
     double eps,
     double h,
     odesolverstate* state,
     ae_state *_state)
{

    _odesolverstate_clear(state);

    ae_assert(n>=1, "ODESolverRKCK: N<1!", _state);
    ae_assert(m>=1, "ODESolverRKCK: M<1!", _state);
    ae_assert(y->cnt>=n, "ODESolverRKCK: Length(Y)<N!", _state);
    ae_assert(x->cnt>=m, "ODESolverRKCK: Length(X)<M!", _state);
    ae_assert(isfinitevector(y, n, _state), "ODESolverRKCK: Y contains infinite or NaN values!", _state);
    ae_assert(isfinitevector(x, m, _state), "ODESolverRKCK: Y contains infinite or NaN values!", _state);
    ae_assert(ae_isfinite(eps, _state), "ODESolverRKCK: Eps is not finite!", _state);
    ae_assert(ae_fp_neq(eps,(double)(0)), "ODESolverRKCK: Eps is zero!", _state);
    ae_assert(ae_isfinite(h, _state), "ODESolverRKCK: H is not finite!", _state);
    odesolver_odesolverinit(0, y, n, x, m, eps, h, state, _state);
}


/*************************************************************************

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool odesolveriteration(odesolverstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double xc;
    double v;
    double h;
    double h2;
    ae_bool gridpoint;
    double err;
    double maxgrowpow;
    ae_int_t klimit;
    ae_bool result;


    
    /*
     * Reverse communication preparations
     * I know it looks ugly, but it works the same way
     * anywhere from C++ to Python.
     *
     * This code initializes locals by:
     * * random values determined during code
     *   generation - on first subroutine call
     * * values from previous call - on subsequent calls
     */
    if( state->rstate.stage>=0 )
    {
        n = state->rstate.ia.ptr.p_int[0];
        m = state->rstate.ia.ptr.p_int[1];
        i = state->rstate.ia.ptr.p_int[2];
        j = state->rstate.ia.ptr.p_int[3];
        k = state->rstate.ia.ptr.p_int[4];
        klimit = state->rstate.ia.ptr.p_int[5];
        gridpoint = state->rstate.ba.ptr.p_bool[0];
        xc = state->rstate.ra.ptr.p_double[0];
        v = state->rstate.ra.ptr.p_double[1];
        h = state->rstate.ra.ptr.p_double[2];
        h2 = state->rstate.ra.ptr.p_double[3];
        err = state->rstate.ra.ptr.p_double[4];
        maxgrowpow = state->rstate.ra.ptr.p_double[5];
    }
    else
    {
        n = 359;
        m = -58;
        i = -919;
        j = -909;
        k = 81;
        klimit = 255;
        gridpoint = ae_false;
        xc = -788;
        v = 809;
        h = 205;
        h2 = -838;
        err = 939;
        maxgrowpow = -526;
    }
    if( state->rstate.stage==0 )
    {
        goto lbl_0;
    }
    
    /*
     * Routine body
     */
    
    /*
     * prepare
     */
    if( state->repterminationtype!=0 )
    {
        result = ae_false;
        return result;
    }
    n = state->n;
    m = state->m;
    h = state->h;
    maxgrowpow = ae_pow(odesolver_odesolvermaxgrow, (double)(5), _state);
    state->repnfev = 0;
    
    /*
     * some preliminary checks for internal errors
     * after this we assume that H>0 and M>1
     */
    ae_assert(ae_fp_greater(state->h,(double)(0)), "ODESolver: internal error", _state);
    ae_assert(m>1, "ODESolverIteration: internal error", _state);
    
    /*
     * choose solver
     */
    if( state->solvertype!=0 )
    {
        goto lbl_1;
    }
    
    /*
     * Cask-Karp solver
     * Prepare coefficients table.
     * Check it for errors
     */
    ae_vector_set_length(&state->rka, 6, _state);
    state->rka.ptr.p_double[0] = (double)(0);
    state->rka.ptr.p_double[1] = (double)1/(double)5;
    state->rka.ptr.p_double[2] = (double)3/(double)10;
    state->rka.ptr.p_double[3] = (double)3/(double)5;
    state->rka.ptr.p_double[4] = (double)(1);
    state->rka.ptr.p_double[5] = (double)7/(double)8;
    ae_matrix_set_length(&state->rkb, 6, 5, _state);
    state->rkb.ptr.pp_double[1][0] = (double)1/(double)5;
    state->rkb.ptr.pp_double[2][0] = (double)3/(double)40;
    state->rkb.ptr.pp_double[2][1] = (double)9/(double)40;
    state->rkb.ptr.pp_double[3][0] = (double)3/(double)10;
    state->rkb.ptr.pp_double[3][1] = -(double)9/(double)10;
    state->rkb.ptr.pp_double[3][2] = (double)6/(double)5;
    state->rkb.ptr.pp_double[4][0] = -(double)11/(double)54;
    state->rkb.ptr.pp_double[4][1] = (double)5/(double)2;
    state->rkb.ptr.pp_double[4][2] = -(double)70/(double)27;
    state->rkb.ptr.pp_double[4][3] = (double)35/(double)27;
    state->rkb.ptr.pp_double[5][0] = (double)1631/(double)55296;
    state->rkb.ptr.pp_double[5][1] = (double)175/(double)512;
    state->rkb.ptr.pp_double[5][2] = (double)575/(double)13824;
    state->rkb.ptr.pp_double[5][3] = (double)44275/(double)110592;
    state->rkb.ptr.pp_double[5][4] = (double)253/(double)4096;
    ae_vector_set_length(&state->rkc, 6, _state);
    state->rkc.ptr.p_double[0] = (double)37/(double)378;
    state->rkc.ptr.p_double[1] = (double)(0);
    state->rkc.ptr.p_double[2] = (double)250/(double)621;
    state->rkc.ptr.p_double[3] = (double)125/(double)594;
    state->rkc.ptr.p_double[4] = (double)(0);
    state->rkc.ptr.p_double[5] = (double)512/(double)1771;
    ae_vector_set_length(&state->rkcs, 6, _state);
    state->rkcs.ptr.p_double[0] = (double)2825/(double)27648;
    state->rkcs.ptr.p_double[1] = (double)(0);
    state->rkcs.ptr.p_double[2] = (double)18575/(double)48384;
    state->rkcs.ptr.p_double[3] = (double)13525/(double)55296;
    state->rkcs.ptr.p_double[4] = (double)277/(double)14336;
    state->rkcs.ptr.p_double[5] = (double)1/(double)4;
    ae_matrix_set_length(&state->rkk, 6, n, _state);
    
    /*
     * Main cycle consists of two iterations:
     * * outer where we travel from X[i-1] to X[i]
     * * inner where we travel inside [X[i-1],X[i]]
     */
    ae_matrix_set_length(&state->ytbl, m, n, _state);
    ae_vector_set_length(&state->escale, n, _state);
    ae_vector_set_length(&state->yn, n, _state);
    ae_vector_set_length(&state->yns, n, _state);
    xc = state->xg.ptr.p_double[0];
    ae_v_move(&state->ytbl.ptr.pp_double[0][0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(j=0; j<=n-1; j++)
    {
        state->escale.ptr.p_double[j] = (double)(0);
    }
    i = 1;
lbl_3:
    if( i>m-1 )
    {
        goto lbl_5;
    }
    
    /*
     * begin inner iteration
     */
lbl_6:
    if( ae_false )
    {
        goto lbl_7;
    }
    
    /*
     * truncate step if needed (beyond right boundary).
     * determine should we store X or not
     */
    if( ae_fp_greater_eq(xc+h,state->xg.ptr.p_double[i]) )
    {
        h = state->xg.ptr.p_double[i]-xc;
        gridpoint = ae_true;
    }
    else
    {
        gridpoint = ae_false;
    }
    
    /*
     * Update error scale maximums
     *
     * These maximums are initialized by zeros,
     * then updated every iterations.
     */
    for(j=0; j<=n-1; j++)
    {
        state->escale.ptr.p_double[j] = ae_maxreal(state->escale.ptr.p_double[j], ae_fabs(state->yc.ptr.p_double[j], _state), _state);
    }
    
    /*
     * make one step:
     * 1. calculate all info needed to do step
     * 2. update errors scale maximums using values/derivatives
     *    obtained during (1)
     *
     * Take into account that we use scaling of X to reduce task
     * to the form where x[0] < x[1] < ... < x[n-1]. So X is
     * replaced by x=xscale*t, and dy/dx=f(y,x) is replaced
     * by dy/dt=xscale*f(y,xscale*t).
     */
    ae_v_move(&state->yn.ptr.p_double[0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->yns.ptr.p_double[0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    k = 0;
lbl_8:
    if( k>5 )
    {
        goto lbl_10;
    }
    
    /*
     * prepare data for the next update of YN/YNS
     */
    state->x = state->xscale*(xc+state->rka.ptr.p_double[k]*h);
    ae_v_move(&state->y.ptr.p_double[0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(j=0; j<=k-1; j++)
    {
        v = state->rkb.ptr.pp_double[k][j];
        ae_v_addd(&state->y.ptr.p_double[0], 1, &state->rkk.ptr.pp_double[j][0], 1, ae_v_len(0,n-1), v);
    }
    state->needdy = ae_true;
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->needdy = ae_false;
    state->repnfev = state->repnfev+1;
    v = h*state->xscale;
    ae_v_moved(&state->rkk.ptr.pp_double[k][0], 1, &state->dy.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    
    /*
     * update YN/YNS
     */
    v = state->rkc.ptr.p_double[k];
    ae_v_addd(&state->yn.ptr.p_double[0], 1, &state->rkk.ptr.pp_double[k][0], 1, ae_v_len(0,n-1), v);
    v = state->rkcs.ptr.p_double[k];
    ae_v_addd(&state->yns.ptr.p_double[0], 1, &state->rkk.ptr.pp_double[k][0], 1, ae_v_len(0,n-1), v);
    k = k+1;
    goto lbl_8;
lbl_10:
    
    /*
     * estimate error
     */
    err = (double)(0);
    for(j=0; j<=n-1; j++)
    {
        if( !state->fraceps )
        {
            
            /*
             * absolute error is estimated
             */
            err = ae_maxreal(err, ae_fabs(state->yn.ptr.p_double[j]-state->yns.ptr.p_double[j], _state), _state);
        }
        else
        {
            
            /*
             * Relative error is estimated
             */
            v = state->escale.ptr.p_double[j];
            if( ae_fp_eq(v,(double)(0)) )
            {
                v = (double)(1);
            }
            err = ae_maxreal(err, ae_fabs(state->yn.ptr.p_double[j]-state->yns.ptr.p_double[j], _state)/v, _state);
        }
    }
    
    /*
     * calculate new step, restart if necessary
     */
    if( ae_fp_less_eq(maxgrowpow*err,state->eps) )
    {
        h2 = odesolver_odesolvermaxgrow*h;
    }
    else
    {
        h2 = h*ae_pow(state->eps/err, 0.2, _state);
    }
    if( ae_fp_less(h2,h/odesolver_odesolvermaxshrink) )
    {
        h2 = h/odesolver_odesolvermaxshrink;
    }
    if( ae_fp_greater(err,state->eps) )
    {
        h = ae_minreal(h2, odesolver_odesolverguaranteeddecay*h, _state);
        goto lbl_6;
    }
    
    /*
     * advance position
     */
    xc = xc+h;
    ae_v_move(&state->yc.ptr.p_double[0], 1, &state->yn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * update H
     */
    h = h2;
    
    /*
     * break on grid point
     */
    if( gridpoint )
    {
        goto lbl_7;
    }
    goto lbl_6;
lbl_7:
    
    /*
     * save result
     */
    ae_v_move(&state->ytbl.ptr.pp_double[i][0], 1, &state->yc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    i = i+1;
    goto lbl_3;
lbl_5:
    state->repterminationtype = 1;
    result = ae_false;
    return result;
lbl_1:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = m;
    state->rstate.ia.ptr.p_int[2] = i;
    state->rstate.ia.ptr.p_int[3] = j;
    state->rstate.ia.ptr.p_int[4] = k;
    state->rstate.ia.ptr.p_int[5] = klimit;
    state->rstate.ba.ptr.p_bool[0] = gridpoint;
    state->rstate.ra.ptr.p_double[0] = xc;
    state->rstate.ra.ptr.p_double[1] = v;
    state->rstate.ra.ptr.p_double[2] = h;
    state->rstate.ra.ptr.p_double[3] = h2;
    state->rstate.ra.ptr.p_double[4] = err;
    state->rstate.ra.ptr.p_double[5] = maxgrowpow;
    return result;
}


/*************************************************************************
ODE solver results

Called after OdeSolverIteration returned False.

INPUT PARAMETERS:
    State   -   algorithm state (used by OdeSolverIteration).

OUTPUT PARAMETERS:
    M       -   number of tabulated values, M>=1
    XTbl    -   array[0..M-1], values of X
    YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]
    Rep     -   solver report:
                * Rep.TerminationType completetion code:
                    * -2    X is not ordered  by  ascending/descending  or
                            there are non-distinct X[],  i.e.  X[i]=X[i+1]
                    * -1    incorrect parameters were specified
                    *  1    task has been solved
                * Rep.NFEV contains number of function calculations

  -- ALGLIB --
     Copyright 01.09.2009 by Bochkanov Sergey
*************************************************************************/
void odesolverresults(odesolverstate* state,
     ae_int_t* m,
     /* Real    */ ae_vector* xtbl,
     /* Real    */ ae_matrix* ytbl,
     odesolverreport* rep,
     ae_state *_state)
{
    double v;
    ae_int_t i;

    *m = 0;
    ae_vector_clear(xtbl);
    ae_matrix_clear(ytbl);
    _odesolverreport_clear(rep);

    rep->terminationtype = state->repterminationtype;
    if( rep->terminationtype>0 )
    {
        *m = state->m;
        rep->nfev = state->repnfev;
        ae_vector_set_length(xtbl, state->m, _state);
        v = state->xscale;
        ae_v_moved(&xtbl->ptr.p_double[0], 1, &state->xg.ptr.p_double[0], 1, ae_v_len(0,state->m-1), v);
        ae_matrix_set_length(ytbl, state->m, state->n, _state);
        for(i=0; i<=state->m-1; i++)
        {
            ae_v_move(&ytbl->ptr.pp_double[i][0], 1, &state->ytbl.ptr.pp_double[i][0], 1, ae_v_len(0,state->n-1));
        }
    }
    else
    {
        rep->nfev = 0;
    }
}


/*************************************************************************
Internal initialization subroutine
*************************************************************************/
static void odesolver_odesolverinit(ae_int_t solvertype,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_int_t m,
     double eps,
     double h,
     odesolverstate* state,
     ae_state *_state)
{
    ae_int_t i;
    double v;

    _odesolverstate_clear(state);

    
    /*
     * Prepare RComm
     */
    ae_vector_set_length(&state->rstate.ia, 5+1, _state);
    ae_vector_set_length(&state->rstate.ba, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 5+1, _state);
    state->rstate.stage = -1;
    state->needdy = ae_false;
    
    /*
     * check parameters.
     */
    if( (n<=0||m<1)||ae_fp_eq(eps,(double)(0)) )
    {
        state->repterminationtype = -1;
        return;
    }
    if( ae_fp_less(h,(double)(0)) )
    {
        h = -h;
    }
    
    /*
     * quick exit if necessary.
     * after this block we assume that M>1
     */
    if( m==1 )
    {
        state->repnfev = 0;
        state->repterminationtype = 1;
        ae_matrix_set_length(&state->ytbl, 1, n, _state);
        ae_v_move(&state->ytbl.ptr.pp_double[0][0], 1, &y->ptr.p_double[0], 1, ae_v_len(0,n-1));
        ae_vector_set_length(&state->xg, m, _state);
        ae_v_move(&state->xg.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,m-1));
        return;
    }
    
    /*
     * check again: correct order of X[]
     */
    if( ae_fp_eq(x->ptr.p_double[1],x->ptr.p_double[0]) )
    {
        state->repterminationtype = -2;
        return;
    }
    for(i=1; i<=m-1; i++)
    {
        if( (ae_fp_greater(x->ptr.p_double[1],x->ptr.p_double[0])&&ae_fp_less_eq(x->ptr.p_double[i],x->ptr.p_double[i-1]))||(ae_fp_less(x->ptr.p_double[1],x->ptr.p_double[0])&&ae_fp_greater_eq(x->ptr.p_double[i],x->ptr.p_double[i-1])) )
        {
            state->repterminationtype = -2;
            return;
        }
    }
    
    /*
     * auto-select H if necessary
     */
    if( ae_fp_eq(h,(double)(0)) )
    {
        v = ae_fabs(x->ptr.p_double[1]-x->ptr.p_double[0], _state);
        for(i=2; i<=m-1; i++)
        {
            v = ae_minreal(v, ae_fabs(x->ptr.p_double[i]-x->ptr.p_double[i-1], _state), _state);
        }
        h = 0.001*v;
    }
    
    /*
     * store parameters
     */
    state->n = n;
    state->m = m;
    state->h = h;
    state->eps = ae_fabs(eps, _state);
    state->fraceps = ae_fp_less(eps,(double)(0));
    ae_vector_set_length(&state->xg, m, _state);
    ae_v_move(&state->xg.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,m-1));
    if( ae_fp_greater(x->ptr.p_double[1],x->ptr.p_double[0]) )
    {
        state->xscale = (double)(1);
    }
    else
    {
        state->xscale = (double)(-1);
        ae_v_muld(&state->xg.ptr.p_double[0], 1, ae_v_len(0,m-1), -1);
    }
    ae_vector_set_length(&state->yc, n, _state);
    ae_v_move(&state->yc.ptr.p_double[0], 1, &y->ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->solvertype = solvertype;
    state->repterminationtype = 0;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&state->y, n, _state);
    ae_vector_set_length(&state->dy, n, _state);
}


void _odesolverstate_init(void* _p, ae_state *_state, ae_bool make_automatic)
{
    odesolverstate *p = (odesolverstate*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_init(&p->yc, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->escale, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->xg, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->dy, 0, DT_REAL, _state, make_automatic);
    ae_matrix_init(&p->ytbl, 0, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->yn, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->yns, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->rka, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->rkc, 0, DT_REAL, _state, make_automatic);
    ae_vector_init(&p->rkcs, 0, DT_REAL, _state, make_automatic);
    ae_matrix_init(&p->rkb, 0, 0, DT_REAL, _state, make_automatic);
    ae_matrix_init(&p->rkk, 0, 0, DT_REAL, _state, make_automatic);
    _rcommstate_init(&p->rstate, _state, make_automatic);
}


void _odesolverstate_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic)
{
    odesolverstate *dst = (odesolverstate*)_dst;
    odesolverstate *src = (odesolverstate*)_src;
    dst->n = src->n;
    dst->m = src->m;
    dst->xscale = src->xscale;
    dst->h = src->h;
    dst->eps = src->eps;
    dst->fraceps = src->fraceps;
    ae_vector_init_copy(&dst->yc, &src->yc, _state, make_automatic);
    ae_vector_init_copy(&dst->escale, &src->escale, _state, make_automatic);
    ae_vector_init_copy(&dst->xg, &src->xg, _state, make_automatic);
    dst->solvertype = src->solvertype;
    dst->needdy = src->needdy;
    dst->x = src->x;
    ae_vector_init_copy(&dst->y, &src->y, _state, make_automatic);
    ae_vector_init_copy(&dst->dy, &src->dy, _state, make_automatic);
    ae_matrix_init_copy(&dst->ytbl, &src->ytbl, _state, make_automatic);
    dst->repterminationtype = src->repterminationtype;
    dst->repnfev = src->repnfev;
    ae_vector_init_copy(&dst->yn, &src->yn, _state, make_automatic);
    ae_vector_init_copy(&dst->yns, &src->yns, _state, make_automatic);
    ae_vector_init_copy(&dst->rka, &src->rka, _state, make_automatic);
    ae_vector_init_copy(&dst->rkc, &src->rkc, _state, make_automatic);
    ae_vector_init_copy(&dst->rkcs, &src->rkcs, _state, make_automatic);
    ae_matrix_init_copy(&dst->rkb, &src->rkb, _state, make_automatic);
    ae_matrix_init_copy(&dst->rkk, &src->rkk, _state, make_automatic);
    _rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic);
}


void _odesolverstate_clear(void* _p)
{
    odesolverstate *p = (odesolverstate*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_clear(&p->yc);
    ae_vector_clear(&p->escale);
    ae_vector_clear(&p->xg);
    ae_vector_clear(&p->y);
    ae_vector_clear(&p->dy);
    ae_matrix_clear(&p->ytbl);
    ae_vector_clear(&p->yn);
    ae_vector_clear(&p->yns);
    ae_vector_clear(&p->rka);
    ae_vector_clear(&p->rkc);
    ae_vector_clear(&p->rkcs);
    ae_matrix_clear(&p->rkb);
    ae_matrix_clear(&p->rkk);
    _rcommstate_clear(&p->rstate);
}


void _odesolverstate_destroy(void* _p)
{
    odesolverstate *p = (odesolverstate*)_p;
    ae_touch_ptr((void*)p);
    ae_vector_destroy(&p->yc);
    ae_vector_destroy(&p->escale);
    ae_vector_destroy(&p->xg);
    ae_vector_destroy(&p->y);
    ae_vector_destroy(&p->dy);
    ae_matrix_destroy(&p->ytbl);
    ae_vector_destroy(&p->yn);
    ae_vector_destroy(&p->yns);
    ae_vector_destroy(&p->rka);
    ae_vector_destroy(&p->rkc);
    ae_vector_destroy(&p->rkcs);
    ae_matrix_destroy(&p->rkb);
    ae_matrix_destroy(&p->rkk);
    _rcommstate_destroy(&p->rstate);
}


void _odesolverreport_init(void* _p, ae_state *_state, ae_bool make_automatic)
{
    odesolverreport *p = (odesolverreport*)_p;
    ae_touch_ptr((void*)p);
}


void _odesolverreport_init_copy(void* _dst, void* _src, ae_state *_state, ae_bool make_automatic)
{
    odesolverreport *dst = (odesolverreport*)_dst;
    odesolverreport *src = (odesolverreport*)_src;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
}


void _odesolverreport_clear(void* _p)
{
    odesolverreport *p = (odesolverreport*)_p;
    ae_touch_ptr((void*)p);
}


void _odesolverreport_destroy(void* _p)
{
    odesolverreport *p = (odesolverreport*)_p;
    ae_touch_ptr((void*)p);
}


#endif

}

