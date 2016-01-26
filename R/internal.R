.onAttach <- 
  function(libname, pkgname) {
    packageStartupMessage(" ****************************************************************************")
    packageStartupMessage(" GRMeta: processing chemical data by David Enot <david.enot@gustaveroussy.fr>")
    packageStartupMessage(" Last Source code available at https://github.com/tonedivad/GRMeta")
    packageStartupMessage(" ****************************************************************************")
  }

.cvf <- function(x, n = 0) ifelse(sum(!is.na(x)) >= n, 100 * 
                                    sd(x, na.rm = T)/mean(x, na.rm = T), NA)
.cvrf <- function(x, n = 0) ifelse(sum(!is.na(x)) >= n, 100 * 
                                     mad(x, na.rm = T)/median(x, na.rm = T), NA)
