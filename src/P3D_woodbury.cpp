#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

bool woodbury_update(const arma::mat & X,
                     const arma::uvec & na,
                     const arma::uvec & nona,
                     arma::mat & out) {

  arma::mat A = X.submat(nona, nona);

  if (na.n_elem == 0) {
    out = A;
    return true;
  }

  arma::mat B = X.submat(nona, na);
  arma::mat C = X.submat(na, nona);
  arma::mat D = X.submat(na, na);

  arma::mat Dinv;
  if (!arma::inv_sympd(Dinv, D)) {
    out.reset();
    return false;
  }

  out = A - (B * Dinv * C);
  return true;
}

}

arma::mat woodbury_cpp(const arma::mat & X,
                       const arma::uvec & na,
                       const arma::uvec & nona) {

  arma::mat out;

  if (!woodbury_update(X, na, nona, out)) {
    Rcpp::stop("woodbury_cpp(): unable to update inverse for dropped rows/columns");
  }

  return out;
}

arma::mat P3D_woodbury_cpp(const arma::mat & Y,
                    const arma::mat & M,
                    const arma::mat & X,
                    const arma::mat & Vinv,
                    const int ncores){

  arma::uword n_marker = M.n_cols;
  arma::vec beta(n_marker, arma::fill::zeros);
  arma::vec betavar(n_marker, arma::fill::zeros);

  #if defined(_OPENMP)
  #pragma omp parallel num_threads(ncores)
  #pragma omp for
  #endif

  for (arma::uword i = 0; i < n_marker; ++i) {

    arma::colvec cur_marker = M.col(i);

    arma::uvec na = find_nonfinite(cur_marker);
    arma::uvec nona = find_finite(cur_marker);

    arma::colvec cur_Y = Y.rows(nona);

    arma::colvec cur_marker_sub = cur_marker.rows(nona);
    arma::mat cur_X = X.rows(nona);
    arma::uword n_fixed = cur_X.n_cols;

    arma::mat cur_Vinv;

    if (Vinv.n_cols == nona.n_elem){
      cur_Vinv = Vinv;
    } else {
      if (!woodbury_update(Vinv, na, nona, cur_Vinv)) {
        beta(i) = arma::datum::nan;
        betavar(i) = arma::datum::nan;
        continue;
      }
    }

    arma::mat XVF = cur_X.t() * cur_Vinv;
    arma::mat fixed_XVX = XVF * cur_X;
    arma::colvec fixed_XVy = XVF * cur_Y;
    arma::colvec Vinv_marker = cur_Vinv * cur_marker_sub;
    arma::colvec marker_XVX_fixed = cur_X.t() * Vinv_marker;
    arma::colvec Vinv_Y = cur_Vinv * cur_Y;

    arma::mat XVX(n_fixed + 1, n_fixed + 1, arma::fill::zeros);
    XVX(0, 0) = arma::dot(cur_marker_sub, Vinv_marker);
    if (n_fixed > 0) {
      XVX.submat(1, 1, n_fixed, n_fixed) = fixed_XVX;
      XVX.submat(1, 0, n_fixed, 0) = marker_XVX_fixed;
      XVX.submat(0, 1, 0, n_fixed) = marker_XVX_fixed.t();
    }

    arma::colvec XVy(n_fixed + 1, arma::fill::zeros);
    XVy(0) = arma::dot(cur_marker_sub, Vinv_Y);
    if (n_fixed > 0) {
      XVy.rows(1, n_fixed) = fixed_XVy;
    }

    arma::mat XVX_i;
    if (!arma::inv_sympd(XVX_i, XVX)) {
      beta(i) = arma::datum::nan;
      betavar(i) = arma::datum::nan;
      continue;
    }

    arma::colvec cur_beta = XVX_i * XVy;

    beta(i) = cur_beta(0);
    betavar(i) = XVX_i(0, 0);

  }

  arma::mat output = join_rows(beta, betavar);
  return output;
}
