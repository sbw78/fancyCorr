#' Function for creating nicely-formatted correlation matrices.
#'
#' \code{fancyCorr} creates nicely-formatted correlation matrices.
#'
#' This function creates correlation matrices. Created as a mercy to my PSYC 5410 students.
#'
#' @param data Data frame containing the variables to be included in the matrix.
#' @param digits The number of decimals to which the correlation coefficeents and
#'  standard deviations (if applicable) should be printed. Defaults to 2.
#' @param p.digits The number of decimals to which the p-values should be printed. Defaults to 3.
#' @param listwise Logical; should all rows containing missing data be removed? Defaults to TRUE.
#' @param p.adj Logical; should the p-values be adjusted for multiple comparisons? Defaults to FALSE.
#' @param sd.on.diag Logical; should the diagonal be replaced with standard deviations? Defaults to TRUE.
#' @param p.on.upper Should p-values for the correlations be printed on the upper triangle? Defaults to TRUE.
#' @param stars Should statistically significant correlation coefficients be starred? Defaults to FALSE.
#' @param star.thresholds A vector of three values corresponding to p-value thresholds that should
#'  receive one, two, or three stars. Only operative if stars=TRUE. Defaults to c(.05, .01, .001).
#' @param method method="pearson" is the default value. The alternatives to be passed to cor are "spearman" and "kendall"
#'
#' @examples
#' set.seed(110)
#' data <- data.frame(matrix(rnorm(50), ncol=5))
#' fancyCorr(data)
#' fancyCorr(data, digits=3, p.digits=3, listwise=FALSE, sd.on.diag=FALSE,
#'   p.adj=FALSE, stars=TRUE, star.thresholds=c(0,0,.1), p.on.upper=FALSE)
#'
#' @export

fancyCorr <- function(data, digits=2, p.digits=3, listwise=TRUE, p.adj=FALSE,
                      sd.on.diag=TRUE, p.on.upper=TRUE, stars=FALSE,
                      star.thresholds=c(.05, .01, .001), method = "pearson") {

  approx.p.threshold <- 1/(10^p.digits)

  # if listwise=T, drop rows containing missing data
  if (listwise==T) {
    data <- data[stats::complete.cases(data),]
  }

  # calculate SDs
  sds <-  apply(data, 2, stats::sd, na.rm=TRUE)

  # calculate correlations, p-values, etc
  corr.matrix <- psych::corr.test(data, method = method)

  # extract and round the correlation coefficients
  corrs <- format(round(corr.matrix$r, digits=digits), nsmall=digits)

  # extract the p-values
  p <- corr.matrix$p

  if (p.adj==T) {
    # replace lower triangle (unadjusted p values) with
    #  the adjusted p values from the upper triangle
    p[lower.tri(p)] <- t(p)[lower.tri(t(p))]
  } else {
    # replace upper triangle (adjusted p values) with
    #  the adjusted p values from the lower triangle
    p[upper.tri(p)] <- t(p)[upper.tri(t(p))]
  }

  # add stars to correlation values if stars=T
  if (stars==T) {
    corrs[p<star.thresholds[1]] <- paste0(corrs[p<star.thresholds[1]], "*")
    corrs[p<star.thresholds[2]] <- paste0(corrs[p<star.thresholds[2]], "*")
    corrs[p<star.thresholds[3]] <- paste0(corrs[p<star.thresholds[3]], "*")
  }

  if (sd.on.diag==T) {
    # replace diagonal with standard deviations in parentheses
    diag(corrs) <- paste0("(", format(sds,digits=digits), ")")
  } else {
    diag(corrs) <- rep("1", times=nrow(corrs))
  }

  if (p.on.upper==T) {
     corrs[upper.tri(corrs)] <-
      ifelse(
        p[upper.tri(p)] < approx.p.threshold,
        paste0("<", approx.p.threshold),
        format(round(p[upper.tri(p)], p.digits), nsmall=p.digits))
  } else {corrs[upper.tri(corrs)] <- ""}

  return(corrs)
}
