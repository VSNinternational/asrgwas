#include <RcppArmadillo.h>
#include <omp.h>

arma::mat schulz_cpp(const arma::mat & X,
                     const arma::mat & Xinv_init,
                     const arma::uvec & na,
                     const int niter) {

  arma::mat Xinv_tmp = Xinv_init;

  arma::mat I = arma::mat(X.n_cols, X.n_cols, arma::fill::eye) * 2;

  Xinv_tmp.rows(na).fill(0);
  Xinv_tmp.cols(na).fill(0);

  for (int i = 1; i <= niter; i++) {

    Xinv_tmp = Xinv_tmp * (I - X * Xinv_tmp);

  }

  return Xinv_tmp;

}

arma::mat P3D_schulz_cpp(const arma::mat & Y,
                        const arma::mat & M,
                        const arma::mat & X,
                        const arma::mat & V,
                        const arma::mat & Vinv,
                        const int niter,
                        const int ncores){

  int n_marker = M.n_cols;
  arma::vec beta(n_marker, arma::fill::zeros);
  arma::vec betavar(n_marker, arma::fill::zeros);

  #if defined(_OPENMP)
  #pragma omp parallel num_threads(ncores)
  #pragma omp for
  #endif

  for (int i = 0; i < n_marker; ++i) {

    arma::mat cur_Y = Y.col(i);

    arma::mat cur_marker = M.col(i);

    arma::uvec na = find_nonfinite(cur_marker);

    arma::mat cur_X = join_rows(cur_marker, X);

    cur_X.rows(na).fill(0);

    arma::mat cur_Vinv;

    if (niter > 0){
      cur_Vinv = schulz_cpp(V, Vinv, na, niter);
    } else{
      cur_Vinv = Vinv;
    }

    arma::mat XV = cur_X.t() * cur_Vinv;

    arma::mat XVX = XV * cur_X;

    arma::mat XVy = XV * cur_Y;

    arma::mat XVX_i;
    arma::inv_sympd(XVX_i, XVX, arma::inv_opts::allow_approx);

    arma::colvec cur_beta = XVX_i * XVy;

    beta(i) = cur_beta(0, 0);

    arma::colvec cur_betavar = arma::diagvec(XVX_i);

    betavar(i) = cur_betavar(0, 0);

  }

  arma::mat output = join_rows(beta, betavar);
  return output;
}
