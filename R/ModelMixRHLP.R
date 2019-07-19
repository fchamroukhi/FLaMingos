#' A Reference Class which represents a fitted MixRHLP model.
#'
#' ModelMixRHLP represents an estimated MixRHLP model.
#'
#' @field param A [ParamMixRHLP][ParamMixRHLP] object. It contains the
#'   estimated values of the parameters.
#' @field stat A [StatMixRHLP][StatMixRHLP] object. It contains all the
#'   statistics associated to the MixRHLP model.
#' @seealso [ParamMixRHLP], [StatMixRHLP]
#' @export
ModelMixRHLP <- setRefClass(
  "ModelMixRHLP",
  fields = list(
    param = "ParamMixRHLP",
    stat = "StatMixRHLP"
  ),
  methods = list(
    plot = function(what = c("estimatedsignal", "regressors", "loglikelihood"), ...) {
      "Plot method.
      \\describe{
        \\item{\\code{what}}{The type of graph requested:
          \\itemize{
            \\item \\code{\"estimatedsignal\" = } Estimated signal (field
              \\code{Ex} of class \\link{StatMixRHLP}).
            \\item \\code{\"regressors\" = } Polynomial regression components
              (fields \\code{polynomials} and \\code{pi_jgk} of class
              \\link{StatMixRHLP}).
            \\item \\code{\"loglikelihood\" = } Value of the log-likelihood for
              each iteration (field \\code{stored_loglik} of class
              \\link{StatMixRHLP}).
          }
        }
        \\item{\\code{\\dots}}{Other graphics parameters.}
      }
      By default, all the above graphs are produced."

      what <- match.arg(what, several.ok = TRUE)

      oldpar <- par(no.readonly = TRUE)
      on.exit(par(oldpar), add = TRUE)

      # yaxislim <- c(min(modelMixRHLP$Y) - 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))), max(modelMixRHLP$Y) + 2 * mean(sqrt(apply(modelMixRHLP$Y, 1, var))))

      colorsvector = rainbow(param$G)

      if (any(what == "estimatedsignal")) {
        # Cluster and means
        nonemptyclusters = length(unique(stat$klas))
        par(mfrow = c(ceiling(sqrt(nonemptyclusters + 1)), round(sqrt(nonemptyclusters + 1))), mai = c(0.6, 0.6, 0.5, 0.25), mgp = c(2, 1, 0))

        matplot(param$fData$X, t(param$fData$Y), type = "l", lty = "solid", col = "black", xlab = "x", ylab = "y", ...)
        title(main = "Dataset")

        for (g in 1:param$G) {
          cluster_g = param$fData$Y[stat$klas == g, , drop = FALSE]

          if (length(cluster_g) != 0) {
            matplot(param$fData$X, t(cluster_g), type = "l", lty = "dotted", col = colorsvector[g], xlab = "x", ylab = "y", ...)
            lines(param$fData$X, stat$Ex[, g, drop = FALSE], col = "black", lty = "solid", lwd = 1.5, ...)
            title(main = sprintf("Cluster %1.1i", g))
          }
        }
      }

      if (any(what == "regressors")) {
        par(mfrow = c(2, 1), mai = c(0.6, 0.8, 0.5, 0.5))
        for (g in 1:param$G) {
          cluster_g = param$fData$Y[stat$klas == g, , drop = FALSE]

          if (length(cluster_g) != 0) {
            matplot(param$fData$X, t(cluster_g), type = "l", lty = "dotted", col = colorsvector[g], xlab = "x", ylab = "y", ...)

            # Polynomial regressors
            for (k in 1:param$K) {
              lines(param$fData$X, stat$polynomials[, k, g], col = "black", lty = "dotted", lwd = 1.5, ...)
            }

            lines(param$fData$X, stat$Ex[, g, drop = FALSE], col = "black", lty = "solid", lwd = 1.5, ...)
            title(main = sprintf("Cluster %1.1i", g))

            matplot(param$fData$X, stat$pi_jgk[1:param$fData$m, , g], type = "l", lty = "solid", xlab = "x", ylab = "Logistic proportions", ylim = c(0, 1), ...)
          }
        }
      }

      if (any(what == "loglikelihood")) {
        par(mfrow = c(1, 1))
        plot.default(1:length(stat$stored_loglik), stat$stored_loglik, type = "l", col = "blue", xlab = "EM iteration number", ylab = "Log-likelihood", ...)
        title(main = "Log-likelihood")
      }

    },

    summary = function(digits = getOption("digits")) {
      "Summary method.
      \\describe{
        \\item{\\code{digits}}{The number of significant digits to use when
          printing.}
      }"

      title <- paste("Fitted mixRHLP model")
      txt <- paste(rep("-", min(nchar(title) + 4, getOption("width"))), collapse = "")

      # Title
      cat(txt)
      cat("\n")
      cat(title)
      cat("\n")
      cat(txt)

      cat("\n")
      cat("\n")
      cat(paste0("MixRHLP model with G = ", param$G, ifelse(param$G > 1, " clusters", " cluster"), " and K = ", param$K, ifelse(param$K > 1, " regimes", " regime"), ":"))
      cat("\n")
      cat("\n")

      tab <- data.frame("log-likelihood" = stat$loglik, "nu" = param$nu,
                        "AIC" = stat$AIC, "BIC" = stat$BIC, "ICL" = stat$ICL,
                        row.names = "", check.names = FALSE)
      print(tab, digits = digits)

      cat("\nClustering table (Number of curves in each clusters):\n")
      print(table(stat$klas))

      cat("\nMixing probabilities (cluster weights):\n")
      pro <- data.frame(param$alpha)
      colnames(pro) <- 1:param$G
      print(pro, digits = digits, row.names = FALSE)

      cat("\n\n")

      txt <- paste(rep("-", min(nchar(title), getOption("width"))), collapse = "")

      for (g in 1:param$G) {
        cat(txt)
        cat("\nCluster ", g, " (G = ", g, "):\n", sep = "")

        cat("\nRegression coefficients:\n\n")
        if (param$p > 0) {
          row.names = c("1", sapply(1:param$p, function(x) paste0("X^", x)))
        } else {
          row.names = "1"
        }

        betas <- data.frame(matrix(param$beta[, , g], nrow = param$p + 1), row.names = row.names)
        colnames(betas) <- sapply(1:param$K, function(x) paste0("Beta(K = ", x, ")"))
        print(betas, digits = digits)

        cat(paste0(ifelse(param$variance_type == "homoskedastic", "\n", "\nVariances:\n\n")))
        sigma2 <- data.frame(t(param$sigma2[, g]))
        if (param$variance_type == "homoskedastic") {
          colnames(sigma2) <- "Sigma2"
          print(sigma2, digits = digits, row.names = FALSE)
        } else {
          colnames(sigma2) = sapply(1:param$K, function(x) paste0("Sigma2(K = ", x, ")"))
          print(sigma2, digits = digits, row.names = FALSE)
        }
        cat("\n")
      }

    }

  )
)