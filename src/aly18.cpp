// clang-format off
#include <cpp4r.hpp>
#include <armadillo4r.hpp>

using namespace cpp4r;
using namespace arma;

namespace {

// Zero-copy Col<double>/Mat<double> wrapping with argument validation
arma::Col<double> as_arma_col(SEXP x, const char* name) {
  if (TYPEOF(x) != REALSXP) {
    cpp4r::stop("'%s' must be a numeric array", name);
  }
  return arma::Col<double>(REAL(x), static_cast<arma::uword>(Rf_xlength(x)), false, false);
}

arma::Mat<double> as_arma_mat(SEXP x, const char* name, int nrow, int ncol) {
  if (TYPEOF(x) != REALSXP) {
    cpp4r::stop("'%s' must be a numeric matrix/array", name);
  }
  if (Rf_xlength(x) != static_cast<R_xlen_t>(nrow) * ncol) {
    cpp4r::stop("'%s' has incompatible dimensions", name);
  }
  return arma::Mat<double>(REAL(x), nrow, ncol, false, false);
}

// Residuals of the ALY (2018) multilateral resistance systems, ported from
// fcn_mrs_ge_ppml.m (mode 0, x = [inward; outward]) and fcn_pstars_ge_ppml.m
// (mode 1, x = [inward; outward; pstar]). Trade costs are t(i, j) with i the
// exporter and j the importer. The last country is the reference one and its
// inward resistance is normalised to one.
arma::Col<double> geppml_residuals(const arma::Col<double>& x, const arma::Mat<double>& t,
                                   const arma::Col<double>& q, const arma::Col<double>& phi,
                                   const arma::Col<double>& beta, double sigma, int mode) {
  const arma::uword n = q.n_elem;
  arma::Col<double> F(x.n_elem);
  arma::Col<double> y1, y2;
  if (mode == 1) {
    arma::Col<double> pq = x.subvec(2 * n, 3 * n - 1) % q;
    y1 = pq / arma::accu(pq);
    y2 = phi % pq / arma::accu(phi % pq);
  } else {
    y1 = q / arma::accu(q);
    y2 = phi % q / arma::accu(phi % q);
  }
  const arma::Col<double> imr = x.subvec(0, n - 1);
  const arma::Col<double> omr = x.subvec(n, 2 * n - 1);
  const arma::Col<double> inward = imr % (t.t() * (y1 % omr));
  const arma::Col<double> outward = omr % (t * (y2 % imr));
  for (arma::uword i = 0; i < n; i++) {
    F(i) = 1 - inward(i);
    F(n + i) = 1 - outward(i);
  }
  F(n - 1) = 1 - x(n - 1);
  if (mode == 1) {
    for (arma::uword i = 0; i < n; i++) {
      F(2 * n + i) = 1 - y1(i) * omr(i) / (std::pow(x(2 * n + i), 1 - sigma) * beta(i));
    }
  }
  if (x.min() <= 0) {
    F += 1e3;
  }
  return F;
}

// Levenberg-Marquardt replacement for MATLAB's fsolve (central finite
// difference Jacobian, iterates are kept strictly positive).
double geppml_solve(arma::Col<double>& x, const arma::Mat<double>& t,
                    const arma::Col<double>& q, const arma::Col<double>& phi,
                    const arma::Col<double>& beta, double sigma, int mode, double tol,
                    int maxiter, bool trace, const char* label) {
  const arma::uword m = x.n_elem;
  arma::Col<double> F = geppml_residuals(x, t, q, phi, beta, sigma, mode);
  double fn = arma::norm(F, 2);
  double lambda = 1e-3;

  for (int iter = 1; iter <= maxiter; iter++) {
    cpp4r::check_user_interrupt();
    if (trace && (iter == 1 || iter % 10 == 0)) {
      Rprintf("%s - iteration: %d - ||F(x)||: %.10g\n", label, iter, fn);
    }
    if (fn < tol) {
      break;
    }

    arma::Mat<double> J(m, m);
    for (arma::uword k = 0; k < m; k++) {
      const double h = 1e-6 * std::max(std::fabs(x(k)), 1e-6);
      arma::Col<double> xp = x, xm = x;
      xp(k) += h;
      xm(k) -= h;
      J.col(k) = (geppml_residuals(xp, t, q, phi, beta, sigma, mode) -
                  geppml_residuals(xm, t, q, phi, beta, sigma, mode)) /
                 (xp(k) - xm(k));
    }

    const arma::Mat<double> JtJ = J.t() * J;
    const arma::Col<double> g = J.t() * F;
    bool accepted = false;
    for (int tries = 0; tries < 30 && !accepted; tries++) {
      arma::Mat<double> A = JtJ;
      A.diag() += lambda * arma::max(JtJ.diag(), arma::Col<double>(m, arma::fill::ones) * 1e-12);
      arma::Col<double> d;
      if (arma::solve(d, A, -g, arma::solve_opts::no_approx)) {
        arma::Col<double> xn = x + d;
        if (xn.min() > 0) {
          arma::Col<double> Fn = geppml_residuals(xn, t, q, phi, beta, sigma, mode);
          double fnn = arma::norm(Fn, 2);
          if (fnn < fn) {
            x = xn;
            F = Fn;
            fn = fnn;
            lambda = std::max(lambda / 10, 1e-15);
            accepted = true;
            continue;
          }
        }
      }
      lambda *= 10;
    }
    if (!accepted) {
      break;
    }
  }

  if (trace) {
    Rprintf("%s - final ||F(x)||: %.10g\n", label, fn);
  }
  if (fn >= tol) {
    cpp4r::warning("%s: no convergence to the requested tolerance (||F(x)|| = %g).", label, fn);
  }
  return fn;
}

SEXP make_vector(const arma::Col<double>& x) {
  SEXP out = PROTECT(Rf_allocVector(REALSXP, x.n_elem));
  std::copy(x.begin(), x.end(), REAL(out));
  UNPROTECT(1);
  return out;
}

}  // namespace

/* roxygen
@title Solve the ALY (2018) general equilibrium with a non-linear solver
@description Port of ge_ppml_indexes.m. Solves for the baseline, conditional and
  full endowment multilateral resistances (inverse of the indexes, as in the
  MATLAB code) with a Levenberg-Marquardt solver instead of fsolve.
@param t_bln matrix of baseline trade costs (exporter in rows, importer in columns).
@param t_ctrf matrix of counterfactual trade costs, same layout as t_bln.
@param q numeric vector of output by country.
@param phi numeric vector of expenditure to output ratios by country.
@param sigma the elasticity of substitution.
@param tol tolerance for the norm of the residuals.
@param maxiter maximum number of iterations per system.
@param trace a boolean indicating whether convergence information should be printed.
@details The last country is the reference country (its inward resistance is
  normalised to one).
*/
[[cpp4r::register]] cpp4r::sexp nls_geppml_(cpp4r::sexp t_bln, cpp4r::sexp t_ctrf,
                                            cpp4r::sexp q, cpp4r::sexp phi, double sigma,
                                            double tol = 1e-12, int maxiter = 1400,
                                            bool trace = true) {
  const arma::Col<double> q_ = as_arma_col(q.data(), "q");
  const arma::Col<double> phi_ = as_arma_col(phi.data(), "phi");
  const int n = static_cast<int>(q_.n_elem);
  if (n < 2 || phi_.n_elem != q_.n_elem) {
    cpp4r::stop("'q' and 'phi' must have the same length (at least 2)");
  }
  const arma::Mat<double> Tb = as_arma_mat(t_bln.data(), "t_bln", n, n);
  const arma::Mat<double> Tc = as_arma_mat(t_ctrf.data(), "t_ctrf", n, n);
  const arma::Col<double> none(n, arma::fill::zeros);

  // I. baseline
  arma::Col<double> x(2 * n, arma::fill::value(0.1));
  double norm_bln = geppml_solve(x, Tb, q_, phi_, none, sigma, 0, tol, maxiter, trace, "Baseline");
  arma::Col<double> imr_bln = x.head(n), omr_bln = x.tail(n);
  const arma::Col<double> beta = (q_ / arma::accu(q_)) % omr_bln;

  // II. conditional
  x.fill(0.1);
  double norm_cdl = geppml_solve(x, Tc, q_, phi_, none, sigma, 0, tol, maxiter, trace, "Conditional");
  arma::Col<double> imr_cdl = x.head(n), omr_cdl = x.tail(n);

  // III. full endowment
  arma::Col<double> xf(3 * n);
  xf.head(2 * n) = x;
  xf.tail(n).ones();
  double norm_full = geppml_solve(xf, Tc, q_, phi_, beta, sigma, 1, tol, maxiter, trace, "Full endowment");

  const char* names[] = {"imr_bln", "omr_bln", "imr_cdl", "omr_cdl", "imr_full", "omr_full",
                         "pstar", "beta", "norm"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
  SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 9));
  SEXP vals[8];
  arma::Col<double> v[8] = {imr_bln, omr_bln, imr_cdl, omr_cdl,
                            xf.subvec(0, n - 1), xf.subvec(n, 2 * n - 1),
                            xf.subvec(2 * n, 3 * n - 1), beta};
  for (int k = 0; k < 8; k++) {
    vals[k] = PROTECT(make_vector(v[k]));
    SET_VECTOR_ELT(out, k, vals[k]);
  }
  arma::Col<double> norms = {norm_bln, norm_cdl, norm_full};
  SEXP norm_sexp = PROTECT(make_vector(norms));
  SET_VECTOR_ELT(out, 8, norm_sexp);
  for (int k = 0; k < 9; k++) {
    SET_STRING_ELT(out_names, k, Rf_mkChar(names[k]));
  }
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  cpp4r::sexp res(out);
  UNPROTECT(11);
  return res;
}
