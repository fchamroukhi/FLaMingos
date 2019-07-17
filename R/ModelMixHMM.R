#' A Reference Class which represents a fitted MixHMM model.
#'
#' ModelMixHMM represents an estimated MixHMM model.
#'
#' @field param A [ParamMixHMM][ParamMixHMM] object. It contains the
#'   estimated values of the parameters.
#' @field stat A [StatMixHMM][StatMixHMM] object. It contains all the
#'   statistics associated to the MixHMM model.
#' @seealso [ParamMixHMM], [StatMixHMM]
#' @export
ModelMixHMM <- setRefClass(
  "ModelMixHMM",
  fields = list(
    param = "ParamMixHMM",
    stat = "StatMixHMM"
  ),
  methods = list(

    plot = function(...) {
      "Plot method
      \\describe{
        \\item{\\code{\\dots}}{Other graphics parameters.}
      }"

      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)

      # yaxislim <- c(min(param$fData$Y) - 2 * mean(sqrt(apply(param$fData$Y, 1, var))), max(param$fData$Y) + 2 * mean(sqrt(apply(param$fData$Y, 1, var))))

      matplot(t(param$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y(t)", ...)
      title(main = "Original time series")

      colorsvec <- rainbow(param$K)
      matplot(t(param$fData$Y), type = "l", lty = "dotted", col = colorsvec[stat$klas], xlab = "x", ylab = "y(t)", ...)
      title(main = "Clustered time series")

      for (k in 1:param$K) {
        if (sum(stat$klas == k) >= 1) {# At least one curve belongs to cluster k

          if (sum(stat$klas == k) == 1) {# Only one curve in cluster k
            matplot(param$fData$Y[stat$klas == k,], type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)", ...)
          } else {
            matplot(t(param$fData$Y[stat$klas == k,]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)", ...)
          }
          title(main = sprintf("Cluster %1.1i", k))
          lines(stat$smoothed[, k], lwd = 1.5, ...)
        }
      }

      plot.default(1:length(stat$stored_loglik), stat$stored_loglik, type = "l", col = "blue", xlab = "EM iteration number", ylab = "Log-likelihood", ...)
    },

    summary = function(digits = getOption("digits")) {
      "Summary method.
      \\describe{
        \\item{\\code{digits}}{The number of significant digits to use when
          printing.}
      }"

      title <- paste("Fitted mixHMM model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("MixHMM model with K = ", param$K,ifelse(param$K > 1, " clusters", " cluster"), " and R = ", param$R, ifelse(param$R > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = stat$loglik, "nu" = param$nu,
                        "AIC" = stat$AIC, "BIC" = stat$BIC,
                        row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table (Number of curves in each clusters):\n")
      print(table(stat$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(t(param$alpha))
      colnames(pro) <- 1:param$K
      print(pro, digits = digits, row.names = FALSE)

      cat("\n\n")

      txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")

      for (k in 1:param$K) {
        cat(txt)
        cat("\nCluster ", k, " (K = ", k, "):\n", sep = "")

        cat("\nMeans:\n\n")
        means <- data.frame(t(param$mu[, k]))
        colnames(means) <- sapply(1:param$R, function(x) paste0("R = ", x))
        print(means, digits = digits, row.names = FALSE)

        cat(paste0(ifelse(param$variance_type == "homoskedastic", "\n", "\nVariances:\n\n")))
        sigma2 <- data.frame(t(param$sigma2[, k]))
        if (param$variance_type == "homoskedastic") {
          colnames(sigma2) <- "Sigma2"
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          colnames(sigma2) = sapply(1:param$R, function(x)
            paste0("R = ", x))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }
  )
)
