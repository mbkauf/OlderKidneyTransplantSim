#' Random draw from joint distribution of KDPI and cold ischemia time
#'
#' @param seed A seed number
#' @param n Number of draws needed
#' @param france KDPI distribution from France?
#'
#' @return two numerical vectors
#' @export
#'
#' @examples
#' values <- rand_kdpi_cold(seed=12345, n=20, france=FALSE)
rand_kdpi_cold <- function(seed=11111, n, france=FALSE) {
  set.seed(seed)
  if (france==TRUE) {
    kdpi_shape1 = 1.045481
    kdpi_shape2 = 0.697808
  } else {
    kdpi_shape1 = 1.458402
    kdpi_shape2 = 1.033602
  }
  myCop <- copula::normalCopula(param=0.130092, dim = 2, dispstr = "un")
  myMvd <- copula::mvdc(copula=myCop, margins=c("gamma", "beta"),
                      paramMargins=list(list(shape = 3.949532,
                                             rate = 0.2128366),
                                        list(shape1 = kdpi_shape1,
                                             shape2 = kdpi_shape2)))
  Z <- copula::rMvdc(n, myMvd)
  colnames(Z) <- c("rec_cold_isch_tm", "kdpi")
  Z[,2] <- Z[,2] * 100
  Z[,2] <- round(Z[,2],0)
  return(Z)
}
