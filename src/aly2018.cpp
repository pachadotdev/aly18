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
vec geppml_residuals(const vec& x, const mat& t, const vec& q, const vec& phi,
                     const vec& beta, double sigma, int mode) {
  const uword n = q.n_elem;
  vec F(x.n_elem);
  vec y1, y2;
  if (mode == 1) {
    const vec pq = x.tail(n) % q;
    y1 = pq / accu(pq);
    y2 = phi % pq / accu(phi % pq);
  } else {
    y1 = q / accu(q);
    y2 = phi % q / accu(phi % q);
  }
  const vec imr = x.head(n);
  const vec omr = x.subvec(n, 2 * n - 1);
  F.head(n) = 1 - imr % (t.t() * (y1 % omr));
  F.subvec(n, 2 * n - 1) = 1 - omr % (t * (y2 % imr));
  F(n - 1) = 1 - x(n - 1);
  if (mode == 1) {
    F.tail(n) = 1 - (y1 % omr) / (pow(x.tail(n), 1 - sigma) % beta);
  }
  if (x.min() <= 0) {
    F += 1e3;
  }
  return F;
}

// Levenberg-Marquardt replacement for MATLAB's fsolve (central finite
// difference Jacobian, iterates are kept strictly positive).
double geppml_solve(vec& x, const mat& t, const vec& q, const vec& phi,
                    const vec& beta, double sigma, int mode, double tol,
                    int maxiter, bool trace, const char* label) {
  const uword m = x.n_elem;
  vec F = geppml_residuals(x, t, q, phi, beta, sigma, mode);
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

    mat J(m, m);
    for (uword k = 0; k < m; k++) {
      const double h = 1e-6 * std::max(std::fabs(x(k)), 1e-6);
      vec xp = x, xm = x;
      xp(k) += h;
      xm(k) -= h;
      J.col(k) = (geppml_residuals(xp, t, q, phi, beta, sigma, mode) -
                  geppml_residuals(xm, t, q, phi, beta, sigma, mode)) /
                 (xp(k) - xm(k));
    }

    const mat JtJ = J.t() * J;
    const vec g = J.t() * F;
    bool accepted = false;
    for (int tries = 0; tries < 30 && !accepted; tries++) {
      mat A = JtJ;
      A.diag() += lambda * max(JtJ.diag(), vec(m, fill::value(1e-12)));
      vec d;
      if (solve(d, A, -g, solve_opts::no_approx)) {
        vec xn = x + d;
        if (xn.min() > 0) {
          vec Fn = geppml_residuals(xn, t, q, phi, beta, sigma, mode);
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
[[cpp4r::register]] list nls_geppml_(doubles_matrix<> t_bln,
                                   doubles_matrix<> t_ctrf,
                                   doubles q, doubles phi, double sigma,
                                   double tol = 1e-12, int maxiter = 1400,
                                   bool trace = true) {
  const vec q_ = as_col(q);
  const vec phi_ = as_col(phi);
  const uword n = q_.n_elem;
  if (n < 2 || phi_.n_elem != q_.n_elem) {
    stop("'q' and 'phi' must have the same length (at least 2)");
  }
  const mat Tb = as_mat(t_bln);
  const mat Tc = as_mat(t_ctrf);
  if (Tb.n_rows != n || Tb.n_cols != n || Tc.n_rows != n || Tc.n_cols != n) {
    stop("'t_bln' and 't_ctrf' must be square matrices matching 'q'");
  }
  const vec none(n, fill::zeros);

  // I. baseline
  vec x(2 * n, fill::value(0.1));
  const double norm_bln = geppml_solve(x, Tb, q_, phi_, none, sigma, 0, tol, maxiter, trace, "Baseline");
  const vec imr_bln = x.head(n), omr_bln = x.tail(n);
  const vec beta = (q_ / accu(q_)) % omr_bln;

  // II. conditional
  x.fill(0.1);
  const double norm_cdl = geppml_solve(x, Tc, q_, phi_, none, sigma, 0, tol, maxiter, trace, "Conditional");
  const vec imr_cdl = x.head(n), omr_cdl = x.tail(n);

  // III. full endowment
  vec xf = join_cols(x, vec(n, fill::ones));
  const double norm_full = geppml_solve(xf, Tc, q_, phi_, beta, sigma, 1, tol, maxiter, trace, "Full endowment");
  const vec imr_full = xf.head(n), omr_full = xf.subvec(n, 2 * n - 1);
  const vec pstar = xf.tail(n);
  const vec norms = {norm_bln, norm_cdl, norm_full};

  using namespace cpp4r::literals;
  return writable::list({
      "imr_bln"_nm = as_doubles(imr_bln),
      "omr_bln"_nm = as_doubles(omr_bln),
      "imr_cdl"_nm = as_doubles(imr_cdl),
      "omr_cdl"_nm = as_doubles(omr_cdl),
      "imr_full"_nm = as_doubles(imr_full),
      "omr_full"_nm = as_doubles(omr_full),
      "pstar"_nm = as_doubles(pstar),
      "beta"_nm = as_doubles(beta),
      "norm"_nm = as_doubles(norms)});
}
