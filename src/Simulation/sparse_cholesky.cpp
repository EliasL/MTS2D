#include "sparse_cholesky.h"
#include "finite.h"
#include <Eigen/SparseCholesky>
#include <stdexcept>
#ifdef MTS_HAVE_CHOLMOD
#include <cholmod.h>
#endif

using Matrix = Eigen::SparseMatrix<double>;
using Vector = Eigen::VectorXd;

struct SparseCholesky::Impl {
  SparseBackend backend;
  Eigen::SimplicialLLT<Matrix> eigen;
  int n = 0;
  bool ready = false;
#ifdef MTS_HAVE_CHOLMOD
  cholmod_common common{};
  cholmod_factor *factor = nullptr;
  cholmod_dense *x = nullptr, *y = nullptr, *e = nullptr;
  Matrix lower; // Lazy copy for nonzero forward coordinate maps only.

  static cholmod_sparse view(const Matrix &matrix) {
    cholmod_sparse a{};
    a.nrow = matrix.rows();
    a.ncol = matrix.cols();
    a.nzmax = matrix.nonZeros();
    a.p = const_cast<int *>(matrix.outerIndexPtr());
    a.i = const_cast<int *>(matrix.innerIndexPtr());
    a.x = const_cast<double *>(matrix.valuePtr());
    a.stype = -1;
    a.itype = CHOLMOD_INT;
    a.xtype = CHOLMOD_REAL;
    a.dtype = CHOLMOD_DOUBLE;
    a.sorted = a.packed = 1;
    return a;
  }

  Vector permute(const Vector &v, bool inverse = false) const {
    const auto *p = static_cast<const int *>(factor->Perm);
    Vector result(n);
    for (int i = 0; i < n; ++i) {
      if (inverse) result[p[i]] = v[i];
      else result[i] = v[p[i]];
    }
    return result;
  }

  Vector solveCholmod(int system, const Vector &rhs) {
    cholmod_dense b{};
    b.nrow = b.nzmax = b.d = n;
    b.ncol = 1;
    b.x = const_cast<double *>(rhs.data());
    b.xtype = CHOLMOD_REAL;
    b.dtype = CHOLMOD_DOUBLE;
    const int ok = cholmod_solve2(system, factor, &b, nullptr, &x, nullptr,
                                 &y, &e, &common);
    if (!ok || common.status != CHOLMOD_OK || !x) {
      throw std::runtime_error("CHOLMOD triangular solve failed.");
    }
    return Eigen::Map<const Vector>(static_cast<const double *>(x->x), n);
  }
#endif

  explicit Impl(SparseBackend selected) : backend(selected) {
    if (!SparseCholesky::available(selected)) {
      throw std::invalid_argument("Requested sparse Cholesky backend is unavailable.");
    }
#ifdef MTS_HAVE_CHOLMOD
    if (backend == SparseBackend::Cholmod) {
      if (!cholmod_start(&common)) throw std::runtime_error("CHOLMOD initialization failed.");
      common.supernodal = CHOLMOD_SUPERNODAL;
      common.final_ll = 1;
      common.print = 0;
      // Compare AMD and nested dissection rather than assuming one is best.
      common.nmethods = 2;
      common.method[0].ordering = CHOLMOD_AMD;
      common.method[1].ordering = CHOLMOD_METIS;
    }
#endif
  }

  ~Impl() {
#ifdef MTS_HAVE_CHOLMOD
    if (backend == SparseBackend::Cholmod) {
      cholmod_free_dense(&x, &common);
      cholmod_free_dense(&y, &common);
      cholmod_free_dense(&e, &common);
      cholmod_free_factor(&factor, &common);
      cholmod_finish(&common);
    }
#endif
  }

  void check(const Vector &v) const {
    if (!ready || v.size() != n) {
      throw std::invalid_argument("Sparse Cholesky factor is not ready or vector has wrong size.");
    }
  }
};

bool SparseCholesky::available(SparseBackend backend) {
  if (backend == SparseBackend::Eigen) return true;
#ifdef MTS_HAVE_CHOLMOD
  return backend == SparseBackend::Cholmod;
#else
  return false;
#endif
}

SparseBackend SparseCholesky::defaultBackend() {
#ifdef MTS_DEFAULT_CHOLMOD
  return SparseBackend::Cholmod;
#else
  return SparseBackend::Eigen;
#endif
}

const char *SparseCholesky::name(SparseBackend backend) {
  return backend == SparseBackend::Eigen ? "eigen" : "cholmod";
}

SparseCholesky::SparseCholesky(SparseBackend backend)
    : impl(std::make_unique<Impl>(backend)) {}
SparseCholesky::~SparseCholesky() = default;

void SparseCholesky::analyzePattern(const Matrix &matrix) {
  impl->ready = false;
  if (matrix.rows() == 0 || matrix.rows() != matrix.cols() || !matrix.isCompressed()) {
    throw std::invalid_argument("Sparse Cholesky requires a nonempty compressed square matrix.");
  }
  impl->n = matrix.rows();
  if (impl->backend == SparseBackend::Eigen) {
    impl->eigen.analyzePattern(matrix);
  }
#ifdef MTS_HAVE_CHOLMOD
  else {
    cholmod_free_factor(&impl->factor, &impl->common);
    auto a = Impl::view(matrix);
    impl->factor = cholmod_analyze(&a, &impl->common);
    if (!impl->factor || impl->common.status != CHOLMOD_OK) {
      throw std::runtime_error("CHOLMOD symbolic analysis failed.");
    }
  }
#endif
}

bool SparseCholesky::factorize(const Matrix &matrix, double shift) {
  impl->ready = false;
  if (impl->n == 0 || matrix.rows() != impl->n || matrix.cols() != impl->n ||
      !matrix.isCompressed() || !finiteValue(shift) || shift < 0) {
    throw std::invalid_argument("Invalid sparse Cholesky factorization input.");
  }
  if (impl->backend == SparseBackend::Eigen) {
    impl->eigen.setShift(shift);
    impl->eigen.factorize(matrix);
    impl->ready = impl->eigen.info() == Eigen::Success;
  }
#ifdef MTS_HAVE_CHOLMOD
  else {
    impl->lower.resize(0, 0);
    auto a = Impl::view(matrix);
    double beta[2] = {shift, 0};
    const int ok = cholmod_factorize_p(&a, beta, nullptr, 0, impl->factor, &impl->common);
    if (impl->common.status < CHOLMOD_OK) {
      throw std::runtime_error("CHOLMOD numerical factorization failed.");
    }
    impl->ready = ok && impl->common.status == CHOLMOD_OK &&
                  impl->factor->minor == impl->factor->n;
  }
#endif
  return impl->ready;
}

Vector SparseCholesky::toTransformed(const Vector &physical) const {
  impl->check(physical);
  // Each minimization starts at its recentered origin, so this is normally zero.
  if (physical.isZero(0)) return Vector::Zero(physical.size());
  if (impl->backend == SparseBackend::Eigen) {
    return impl->eigen.matrixU() * (impl->eigen.permutationP() * physical);
  }
#ifdef MTS_HAVE_CHOLMOD
  if (impl->lower.rows() == 0) {
    // Conversion consumes the numerical factor: convert a copy, keeping the
    // original supernodal factor and its reusable solve workspaces intact.
    auto *copy = cholmod_copy_factor(impl->factor, &impl->common);
    if (!copy) throw std::runtime_error("CHOLMOD factor copy failed.");
    auto *l = cholmod_factor_to_sparse(copy, &impl->common);
    cholmod_free_factor(&copy, &impl->common);
    if (!l) throw std::runtime_error("CHOLMOD factor conversion failed.");
    const auto *p = static_cast<const int *>(l->p);
    impl->lower = Eigen::Map<const Matrix>(impl->n, impl->n, p[impl->n], p,
        static_cast<const int *>(l->i), static_cast<const double *>(l->x));
    cholmod_free_sparse(&l, &impl->common);
  }
  return impl->lower.transpose() * impl->permute(physical);
#else
  throw std::logic_error("Unavailable CHOLMOD backend.");
#endif
}

Vector SparseCholesky::toPhysical(const Vector &z) const {
  impl->check(z);
  if (impl->backend == SparseBackend::Eigen) {
    return impl->eigen.permutationPinv() * impl->eigen.matrixU().solve(z);
  }
#ifdef MTS_HAVE_CHOLMOD
  return impl->permute(impl->solveCholmod(CHOLMOD_Lt, z), true);
#else
  throw std::logic_error("Unavailable CHOLMOD backend.");
#endif
}

Vector SparseCholesky::transformGradient(const Vector &g) const {
  impl->check(g);
  if (impl->backend == SparseBackend::Eigen) {
    return impl->eigen.matrixL().solve(impl->eigen.permutationP() * g);
  }
#ifdef MTS_HAVE_CHOLMOD
  return impl->solveCholmod(CHOLMOD_L, impl->permute(g));
#else
  throw std::logic_error("Unavailable CHOLMOD backend.");
#endif
}

Vector SparseCholesky::solve(const Vector &rhs) const {
  impl->check(rhs);
  if (impl->backend == SparseBackend::Eigen) return impl->eigen.solve(rhs);
#ifdef MTS_HAVE_CHOLMOD
  return impl->solveCholmod(CHOLMOD_A, rhs);
#else
  throw std::logic_error("Unavailable CHOLMOD backend.");
#endif
}
