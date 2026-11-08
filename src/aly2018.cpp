// clang-format off
#include <cpp4r.hpp>
#include <armadillo4r.hpp>

using namespace cpp4r;
using namespace arma;

namespace {

// Residuals of the ALY (2018) multilateral resistance systems, ported from
// fcn_mrs_ge_ppml.m (mode 0, x = [inward; outward]) and fcn_pstars_ge_ppml.m
// (mode 1, x = [inward; outward; pstar]). Trade costs are t(i, j) with i the
// exporter and j the importer. The last country is the reference one and its
// inward resistance is normalised to one.
Col<double> geppml_residuals(const Col<double>& x, const Mat<double>& t,
                                   const Col<double>& q, const Col<double>& phi,
                                   const Col<double>& beta, double sigma, int mode) {
  const uword n = q.n_elem;
  Col<double> F(x.n_elem);
  Col<double> y1, y2;
  if (mode == 1) {
    Col<double> pq = x.subvec(2 * n, 3 * n - 1) % q;
    y1 = pq / accu(pq);
    y2 = phi % pq / accu(phi % pq);
  } else {
    y1 = q / accu(q);
    y2 = phi % q / accu(phi % q);
  }
  const Col<double> imr = x.subvec(0, n - 1);
  const Col<double> omr = x.subvec(n, 2 * n - 1);
  const Col<double> inward = imr % (t.t() * (y1 % omr));
  const Col<double> outward = omr % (t * (y2 % imr));
  for (uword i = 0; i < n; i++) {
    F(i) = 1 - inward(i);
    F(n + i) = 1 - outward(i);
  }
  F(n - 1) = 1 - x(n - 1);
  if (mode == 1) {
    for (uword i = 0; i < n; i++) {
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
double geppml_solve(Col<double>& x, const Mat<double>& t,
                    const Col<double>& q, const Col<double>& phi,
                    const Col<double>& beta, double sigma, int mode, double tol,
                    int maxiter, bool trace, const char* label) {
  const uword m = x.n_elem;
  Col<double> F = geppml_residuals(x, t, q, phi, beta, sigma, mode);
  double fn = norm(F, 2);
  double lambda = 1e-3;

  for (int iter = 1; iter <= maxiter; iter++) {
    check_user_interrupt();
    if (trace && (iter == 1 || iter % 10 == 0)) {
      Rprintf("%s - iteration: %d - ||F(x)||: %.10g\n", label, iter, fn);
    }
    if (fn < tol) {
      break;
    }

    Mat<double> J(m, m);
    for (uword k = 0; k < m; k++) {
      const double h = 1e-6 * std::max(std::fabs(x(k)), 1e-6);
      Col<double> xp = x, xm = x;
      xp(k) += h;
      xm(k) -= h;
      J.col(k) = (geppml_residuals(xp, t, q, phi, beta, sigma, mode) -
                  geppml_residuals(xm, t, q, phi, beta, sigma, mode)) /
                 (xp(k) - xm(k));
    }

    const Mat<double> JtJ = J.t() * J;
    const Col<double> g = J.t() * F;
    bool accepted = false;
    for (int tries = 0; tries < 30 && !accepted; tries++) {
      Mat<double> A = JtJ;
      A.diag() += lambda * max(JtJ.diag(), Col<double>(m, fill::ones) * 1e-12);
      Col<double> d;
      if (solve(d, A, -g, solve_opts::no_approx)) {
        Col<double> xn = x + d;
        if (xn.min() > 0) {
          Col<double> Fn = geppml_residuals(xn, t, q, phi, beta, sigma, mode);
          double fnn = norm(Fn, 2);
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
    warning("%s: no convergence to the requested tolerance (||F(x)|| = %g).", label, fn);
  }
  return fn;
}

SEXP make_vector(const Col<double>& x) {
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
[[cpp4r::register]] sexp nls_geppml_(sexp t_bln, sexp t_ctrf,
                                            sexp q, sexp phi, double sigma,
                                            double tol = 1e-12, int maxiter = 1400,
                                            bool trace = true) {
  const Col<double> q_ = as_Col(q.data());
  const Col<double> phi_ = as_Col(phi.data());
  const int n = static_cast<int>(q_.n_elem);
  if (n < 2 || phi_.n_elem != q_.n_elem) {
    stop("'q' and 'phi' must have the same length (at least 2)");
  }
  const Mat<double> Tb = as_Mat(t_bln.data());
  const Mat<double> Tc = as_Mat(t_ctrf.data());
  if (Tb.n_rows != n || Tb.n_cols != n || Tc.n_rows != n || Tc.n_cols != n) {
    stop("'t_bln' and 't_ctrf' must be square matrices matching 'q'");
  }
  const Col<double> none(n, fill::zeros);

  // I. baseline
  Col<double> x(2 * n, fill::value(0.1));
  double norm_bln = geppml_solve(x, Tb, q_, phi_, none, sigma, 0, tol, maxiter, trace, "Baseline");
  Col<double> imr_bln = x.head(n), omr_bln = x.tail(n);
  const Col<double> beta = (q_ / accu(q_)) % omr_bln;

  // II. conditional
  x.fill(0.1);
  double norm_cdl = geppml_solve(x, Tc, q_, phi_, none, sigma, 0, tol, maxiter, trace, "Conditional");
  Col<double> imr_cdl = x.head(n), omr_cdl = x.tail(n);

  // III. full endowment
  Col<double> xf(3 * n);
  xf.head(2 * n) = x;
  xf.tail(n).ones();
  double norm_full = geppml_solve(xf, Tc, q_, phi_, beta, sigma, 1, tol, maxiter, trace, "Full endowment");

  const char* names[] = {"imr_bln", "omr_bln", "imr_cdl", "omr_cdl", "imr_full", "omr_full",
                         "pstar", "beta", "norm"};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 9));
  SEXP out_names = PROTECT(Rf_allocVector(STRSXP, 9));
  SEXP vals[8];
  Col<double> v[8] = {imr_bln, omr_bln, imr_cdl, omr_cdl,
                            xf.subvec(0, n - 1), xf.subvec(n, 2 * n - 1),
                            xf.subvec(2 * n, 3 * n - 1), beta};
  for (int k = 0; k < 8; k++) {
    vals[k] = PROTECT(make_vector(v[k]));
    SET_VECTOR_ELT(out, k, vals[k]);
  }
  Col<double> norms = {norm_bln, norm_cdl, norm_full};
  SEXP norm_sexp = PROTECT(make_vector(norms));
  SET_VECTOR_ELT(out, 8, norm_sexp);
  for (int k = 0; k < 9; k++) {
    SET_STRING_ELT(out_names, k, Rf_mkChar(names[k]));
  }
  Rf_setAttrib(out, R_NamesSymbol, out_names);
  sexp res(out);
  UNPROTECT(11);
  return res;
}
