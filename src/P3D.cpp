#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

arma::mat P3D_cpp(const arma::mat & Y,
                  const arma::mat & M,
                  const arma::mat & X,
                  const arma::mat & Vinv,
                  const int ncores){

  arma::uword n_marker = M.n_cols;
  arma::uword n_fixed = X.n_cols;
  arma::vec beta(n_marker, arma::fill::zeros);
  arma::vec betavar(n_marker, arma::fill::zeros);

  arma::mat Xt = X.t();
  arma::colvec Vinv_Y = Vinv * Y;
  arma::mat XVF = Xt * Vinv;
  arma::mat fixed_XVX = XVF * X;
  arma::colvec fixed_XVy = Xt * Vinv_Y;

  #if defined(_OPENMP)
  #pragma omp parallel num_threads(ncores)
  #pragma omp for
  #endif

  for (arma::uword i = 0; i < n_marker; ++i) {

    arma::colvec cur_marker = M.col(i);
    arma::colvec Vinv_marker = Vinv * cur_marker;
    arma::colvec marker_XVX_fixed = Xt * Vinv_marker;

    arma::mat XVX(n_fixed + 1, n_fixed + 1, arma::fill::zeros);
    XVX(0, 0) = arma::dot(cur_marker, Vinv_marker);
    if (n_fixed > 0) {
      XVX.submat(1, 1, n_fixed, n_fixed) = fixed_XVX;
      XVX.submat(1, 0, n_fixed, 0) = marker_XVX_fixed;
      XVX.submat(0, 1, 0, n_fixed) = marker_XVX_fixed.t();
    }

    arma::colvec XVy(n_fixed + 1, arma::fill::zeros);
    XVy(0) = arma::dot(cur_marker, Vinv_Y);
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
