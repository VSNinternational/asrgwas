P3D_cpp <- function(Y, M, X, Vinv, ncores) {
    .Call(`_ASRgwas_P3D_cpp`, Y, M, X, Vinv, ncores)
}
schulz_cpp <- function(X, Xinv_init, na, niter) {
    .Call(`_ASRgwas_schulz_cpp`, X, Xinv_init, na, niter)
}
P3D_schulz_cpp <- function(Y, M, X, V, Vinv, niter, ncores) {
    .Call(`_ASRgwas_P3D_schulz_cpp`, Y, M, X, V, Vinv, niter, ncores)
}
woodbury_cpp <- function(X, na, nona) {
    .Call(`_ASRgwas_woodbury_cpp`, X, na, nona)
}
P3D_woodbury_cpp <- function(Y, M, X, Vinv, ncores) {
    .Call(`_ASRgwas_P3D_woodbury_cpp`, Y, M, X, Vinv, ncores)
}
