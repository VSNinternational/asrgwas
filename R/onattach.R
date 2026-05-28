.onAttach <- function(libname, pkgname){
  ## This has been moved here so the environment is only set once.
  vsni.init()
  packageStartupMessage(pkgname , " attached from ", libname )
}
vsni.init <- function(askuser = FALSE, msgout = message) {
  ## asreml check
  if (requireNamespace("asreml", quietly = TRUE)) {
    if (packageVersion("asreml") < "4.1.0.176") {
      packageStartupMessage(
        "Incompatible version of asreml detected (>= 4.1.0.176 required)."
      )
    }
  }
  TRUE
}
