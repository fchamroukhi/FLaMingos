#' A Reference Class which represents a fitted MixHMMR model.
#'
#' ModelMixHMMR represents an estimated MixHMMR model.
#'
#' @field param A [ParamMixHMMR][ParamMixHMMR] object. It contains the
#'   estimated values of the parameters.
#' @field stat A [StatMixHMMR][StatMixHMMR] object. It contains all the
#'   statistics associated to the MixHMMR model.
#' @seealso [ParamMixHMMR], [StatMixHMMR]
#' @export
ModelMixHMMR <- setRefClass(
  "ModelMixHMMR",
  fields = list(
    param = "ParamMixHMMR",
    stat = "StatMixHMMR"
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

      matplot(param$fData$X, t(param$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y(t)", ...)
      title(main = "Original time series")

      colorsvec <- rainbow(param$K)
      matplot(param$fData$X, t(param$fData$Y), type = "l", lty = "dotted", col = colorsvec[stat$klas], xlab = "x", ylab = "y(t)", ...)
      title(main = "Clustered time series")

      for (k in 1:param$K) {
        if (sum(stat$klas == k) >= 1) {# At least one curve belongs to cluster k

          if (sum(stat$klas == k) == 1) {# Only one curve in cluster k
            matplot(param$fData$X, param$fData$Y[stat$klas == k,], type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)", ...)
          } else {
            matplot(param$fData$X, t(param$fData$Y[stat$klas == k,]), type = "l", lty = "dotted", col = colorsvec[k], xlab = "x", ylab = "y(t)", ...)
          }
          title(main = sprintf("Cluster %1.1i", k))
          lines(param$fData$X, stat$smoothed[, k], lwd = 1.5, ...)
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

      title <- paste("Fitted mixHMMR model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(
        paste0("MixHMMR model with K = ", param$K, ifelse(param$K > 1, " clusters", " cluster"), " and R = ", param$R, ifelse(param$R > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = stat$loglik, "nu" = param$nu,
          "AIC" = stat$AIC, "BIC" = stat$BIC, "ICL" = stat$ICL1,
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

        cat("\nRegression coefficients:\n\n")
        if (param$p > 0) {
          row.names = c("1", sapply(1:param$p, function(x) paste0("X^", x)))
        } else {
          row.names = "1"
        }

        if (param$p > 0) {
          row.names = c("1", sapply(1:param$p, function(x) paste0("X^", x)))
          betas <- data.frame(param$beta[, , k], row.names = row.names)
        } else {
          row.names = "1"
          betas <- data.frame(t(param$beta[, , k]), row.names = row.names)
        }
        colnames(betas) <- sapply(1:param$R, function(x) paste0("Beta(R = ", x, ")"))
        print(betas, digits = digits)

        cat(paste0(ifelse(param$variance_type == "homoskedastic", "\n", "\nVariances:\n\n")))
        sigma2 <- data.frame(t(param$sigma2[, k]))
        if (param$variance_type == "homoskedastic") {
          colnames(sigma2) <- "Sigma2"
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          colnames(sigma2) = sapply(1:param$R, function(x) paste0("Sigma2(R = ", x, ")"))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }
  )
)
