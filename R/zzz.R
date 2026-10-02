.onAttach <- function(libname, pkgname) {
  pkg_env <- asNamespace(pkgname)
  if (!exists('internal_calibration_maps', envir = pkg_env)) {
    packageStartupMessage(
      '-------------------------------------------------------------------\n',
      'WARNING: \'internal_calibration_maps\' was not found in sysdata.rda.\n',
      'Please run \'data-raw/nhanes_calibration_matrices.R\' using \'nhanesR\'\n',
      'to compile the required NHANES reference matrices before calibration.\n',
      '-------------------------------------------------------------------'
    )
  }
}
