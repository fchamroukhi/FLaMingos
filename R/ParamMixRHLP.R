#' A Reference Class which contains parameters of a mixture of RHLP models.
#'
#' ParamMixRHLP contains all the parameters of a mixture of RHLP models.
#'
#' @field fData [FData][FData] object representing the sample (covariates/inputs
#'   `X` and observed responses/outputs `Y`).
#' @field K The number of clusters (Number of RHLP models).
#' @field R The number of regimes (RHLP components) for each cluster.
#' @field p The order of the polynomial regression.
#' @field q The dimension of the logistic regression. For the purpose of
#'   segmentation, it must be set to 1.
#' @field variance_type Character indicating if the model is homoskedastic
#'   (`variance_type = "homoskedastic"`) or heteroskedastic (`variance_type =
#'   "heteroskedastic"`). By default the model is heteroskedastic.
#' @field alpha Cluster weights. Matrix of dimension \eqn{(1, K)}.
#' @field W Parameters of the logistic process. \eqn{\boldsymbol{W} =
#'   (\boldsymbol{w}_{1},\dots,\boldsymbol{w}_{K})}{W = (w_{1},\dots,w_{K})} is
#'   an array of dimension \eqn{(q + 1, R - 1, K)}, with \eqn{\boldsymbol{w}_{g}
#'   = (\boldsymbol{w}_{g,1},\dots,\boldsymbol{w}_{g,R-1})}{w_{g} =
#'   (w_{g,1},\dots,w_{g,R-1})}, \eqn{g = 1,\dots,K}, and `q` the order of the
#'   logistic regression. `q` is fixed to 1 by default.
#' @field beta Parameters of the polynomial regressions. \eqn{\boldsymbol{\beta}
#'   = (\boldsymbol{\beta}_{1},\dots,\boldsymbol{\beta}_{K})}{\beta =
#'   (\beta_{1},\dots,\beta_{K})} is an array of dimension \eqn{(p + 1, R, K)},
#'   with \eqn{\boldsymbol{\beta}_{g} =
#'   (\boldsymbol{\beta}_{g,1},\dots,\boldsymbol{\beta}_{g,R})}{\beta_{g} =
#'   (\beta_{g,1},\dots,\beta_{g,R})}, \eqn{g = 1,\dots,K}, `p` the order of the
#'   polynomial regression. `p` is fixed to 3 by default.
#' @field sigma2 The variances for the `K` clusters. If MixRHLP model is
#'   heteroskedastic (`variance_type = "heteroskedastic"`) then `sigma2` is a
#'   matrix of size \eqn{(R, K)} (otherwise MixRHLP model is homoskedastic
#'   (`variance_type = "homoskedastic"`) and `sigma2` is a matrix of size
#'   \eqn{(K, 1)}).
#' @field nu The degree of freedom of the MixRHLP model representing the
#'   complexity of the model.
#' @field phi A list giving the regression design matrices for the polynomial
#'   and the logistic regressions.
#' @export
ParamMixRHLP <- setRefClass(
  "ParamMixRHLP",
  fields = list(
    fData = "FData",
    phi = "list",

    K = "numeric", # Number of clusters
    R = "numeric", # Number of regimes
    p = "numeric", # Dimension of beta (order of polynomial regression)
    q = "numeric", # Dimension of w (order of logistic regression)
    variance_type = "character",
    nu = "numeric", # Degree of freedom

    alpha = "matrix", # Cluster weights
    W = "array", # W = (W1,...WK), Wg = (Wg1,...,w_gR-1) parameters of the logistic process: matrix of dimension [(q+1)x(R-1)] with q the order of logistic regression.
    beta = "array", # beta = (beta_1,...,beta_K), beta_g = (beta_g1,...,beta_gR) polynomial regression coefficients: matrix of dimension [(p+1)xR] p being the polynomial degree.
    sigma2 = "matrix" # sigma2 = (sigma_g1,...,sigma_gR) : the variances for the R regmies.
  ),
  methods = list(
    initialize = function(fData = FData(numeric(1), matrix(1)), K = 1, R = 1, p = 3, q = 1, variance_type = "heteroskedastic") {

      fData <<- fData

      phi <<- designmatrix(x = fData$X, p = p, q = q, n = fData$n)

      K <<- K
      R <<- R
      p <<- p
      q <<- q
      variance_type <<- variance_type

      if (variance_type == "homoskedastic") {
        nu <<- (K - 1) + K * ((q + 1) * (R - 1) + R * (p + 1) + 1)
      } else {
        nu <<- (K - 1) + K * ((q + 1) * (R - 1) + R * (p +  1) + R)
      }

      W <<- array(0, dim = c(q + 1, R - 1, K))
      beta <<- array(NA, dim = c(p + 1, R, K))
      if (variance_type == "homoskedastic") {
        sigma2 <<- matrix(NA, K)
      } else {
        sigma2 <<- matrix(NA, R, K)
      }
      alpha <<- matrix(NA, K)
    },

    initParam = function(init_kmeans = TRUE, try_algo = 1) {
      "Method to initialize parameters \\code{alpha}, \\code{W}, \\code{beta}
      and \\code{sigma2}.

      If \\code{init_kmeans = TRUE} then the curve partition is initialized by
      the R-means algorithm. Otherwise the curve partition is initialized
      randomly.

      If \\code{try_algo = 1} then \\code{beta} and \\code{sigma2} are
      initialized by segmenting  the time series \\code{Y} uniformly into
      \\code{R} contiguous segments. Otherwise, \\code{W}, \\code{beta} and
      \\code{sigma2} are initialized by segmenting randomly the time series
      \\code{Y} into \\code{R} segments."

      # 1. Initialization of cluster weights
      alpha <<- 1 / (K * ones(K, 1))

      # 2. Initialization of the model parameters for each cluster: W, betak and sigmak

      # Setting W
      if (try_algo == 1) {
        for (g in (1:K)) { # Random initialization of parameter vector for IRLS
          W[, , g] <<- zeros(q + 1, R - 1)
        }
      } else {
        for (g in (1:K)) { # Random initialization of parameter vector for IRLS
          W[, , g] <<- rand(q + 1, R - 1)
        }
      }

      # betagk and sigma2_gk
      if (init_kmeans) {

        kmeans_res <- kmeans(fData$Y, K, nbr_runs = 20, nbr_iter_max = 400, verbose = FALSE)
        klas <- kmeans_res$klas
        for (g in 1:K) {
          Xg <- fData$Y[klas == g,]
          initRegressionParam(Xg, g, try_algo)
        }
      } else {
        ind <- sample(fData$n)
        for (g in 1:K) {
          if (g < K) {
            Xg <- fData$Y[ind[((g - 1) * round(fData$n / K) + 1):(g * round(fData$n / K))],]
          } else {
            Xg <- fData$Y[ind[((g - 1) * round(fData$n / K) + 1):length(ind)],]
          }
          initRegressionParam(Xg, g, try_algo)
        }
      }
    },

    initRegressionParam = function(Xg, g, try_algo  = 1) {
      "Initialize the matrix of polynomial regression coefficients beta_g for
      the cluster \\code{g}."
      n <- nrow(Xg)
      m <- ncol(Xg)
      if (try_algo == 1) { # Uniform segmentation into R contiguous segments, and then a regression
        zi <- round(m / R) - 1

        beta_k <- matrix(NA, p + 1, R)
        sigma <- c()

        for (k in 1:R) {
          i <- (k - 1) * zi + 1
          j <- k * zi
          Xij <- Xg[, i:j, drop = FALSE]
          Xij <- matrix(t(Xij), ncol = 1)
          phi_ij <- phi$XBeta[i:j, , drop = FALSE]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Xij
          beta_k[, k] <- bk

          if (variance_type == "homoskedastic") {
            sigma <- var(Xij)
          } else {
            mk <- j - i + 1 # length(Xij);
            z <- Xij - Phi_ij %*% bk

            sk <- t(z) %*% z / (n * mk)

            sigma[k] <- sk

          }
        }
      } else {# Random segmentation into R contiguous segments, and then a regression
        Lmin <- round(m / R) # nbr pts min into one segment
        tk_init <- zeros(1, R + 1)
        R_1 <- R
        for (k in 2:R) {
          R_1 <- R_1 - 1

          temp <- tk_init[k - 1] + Lmin:(m - (R_1 * Lmin) - tk_init[k - 1])

          ind <- sample(length(temp))

          tk_init[k] <- temp[ind[1]]
        }
        tk_init[R + 1] <- m
        beta_k <- matrix(NA, p + 1, R)
        sigma <- c()
        for (k in 1:R) {
          i <- tk_init[k] + 1
          j <- tk_init[k + 1]
          Xij <- Xg[, i:j]
          Xij <- matrix(t(Xij), ncol = 1)
          phi_ij <- phi$XBeta[i:j, , drop = FALSE]
          Phi_ij <- repmat(phi_ij, n, 1)

          bk <- solve(t(Phi_ij) %*% Phi_ij) %*% t(Phi_ij) %*% Xij
          beta_k[, k] <- bk

          if (variance_type == "homoskedastic") {
            sigma <- var(Xij)
          } else {
            mk <- j - i + 1 #length(Xij);
            z <- Xij - Phi_ij %*% bk

            sk <- t(z) %*% z / (n * mk)

            sigma[k] <- sk

          }
        }
      }

      beta[, , g] <<- beta_k
      if (variance_type == "homoskedastic") {
        sigma2[g] <<- sigma
      } else {
        sigma2[, g] <<- sigma
      }
    },

    CMStep = function(statMixRHLP, verbose_IRLS = FALSE) {
      "Method which implements the M-step of the CEM algorithm to learn the
      parameters of the MixRHLP model based on statistics provided by the
      object \\code{statMixRHLP} of class \\link{StatMixRHLP} (which contains
      the E-step and the C-step)."

      good_segmentation = TRUE

      alpha <<- t(colSums(statMixRHLP$c_ig)) / fData$n

      # Maximization w.r.t betagk et sigma2_gk
      cluster_labels <- t(repmat(statMixRHLP$klas, 1, fData$m)) # [m x n]
      cluster_labels <- as.vector(cluster_labels)

      for (g in 1:K) {
        Xg = fData$vecY[cluster_labels == g,] # Cluster g (found from a hard clustering)
        tauijk <- as.matrix(statMixRHLP$tau_ijgk[cluster_labels == g, , g]) #[(ng xm) x R]

        if (variance_type == "homoskedastic") {
          s <- 0
        } else {
          sigma_gk <- zeros(R, 1)
        }

        beta_gk <- matrix(NA, p + 1, R)

        for (k in 1:R) {

          segments_weights <- tauijk[, k, drop = F]
          phigk <- (sqrt(segments_weights) %*% ones(1, p + 1)) * phi$XBeta[cluster_labels == g,] # [(ng*m)*(p+1)]
          Xgk <- sqrt(segments_weights) * Xg

          # maximization w.r.t beta_gk: Weighted least squares
          beta_gk[, k] <- solve(t(phigk) %*% phigk + .Machine$double.eps * diag(p + 1)) %*% t(phigk) %*% Xgk # Maximization w.r.t betagk

          #    the same as
          #                 W_gk = diag(cluster_weights.*segment_weights);
          #                 beta_gk(:,k) = inv(phiBeta'*W_gk*phiBeta)*phiBeta'*W_gk*X;
          #   Maximization w.r.t au sigma_gk :
          if (variance_type == "homoskedastic") {
            sk <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2)
            s <- s + sk
            sigma_gk <- s / sum(tauijk)
          } else {
            sigma_gk[k] <-
              colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2) / (sum(segments_weights))
            if ((sum(segments_weights) == 0)) {
              good_segmentation = FALSE
              return(list(0, good_segmentation))
            }
          }
        }

        beta[, , g] <<- beta_gk
        if (variance_type == "homoskedastic") {
          sigma2[g] <<- sigma_gk
        } else {
          sigma2[, g] <<- sigma_gk

        }

        # Maximization w.r.t W

        # Setting of W[,,g]
        Wg_init <- matrix(W[, , g], nrow = q + 1)

        res_irls <- IRLS(phi$Xw[cluster_labels == g,], tauijk, ones(nrow(tauijk), 1), Wg_init, verbose_IRLS)

        W[, , g] <<- res_irls$W
        piik <- res_irls$piik
        reg_irls <- res_irls$reg_irls
      }
      return(list(reg_irls, good_segmentation))
    },

    MStep = function(statMixRHLP, verbose_IRLS = FALSE) {
      "Method which implements the M-step of the EM algorithm to learn the
      parameters of the MixRHLP model based on statistics provided by the
      object \\code{statMixRHLP} of class \\link{StatMixRHLP} (which contains
      the E-step)."

      alpha <<- t(colSums(statMixRHLP$h_ig)) / fData$n
      for (g in 1:K) {
        temp <- repmat(statMixRHLP$h_ig[, g], 1, fData$m) # [m x n]
        cluster_weights <-
          matrix(t(temp), fData$m * fData$n, 1) # Cluster_weights(:) [mn x 1]
        tauijk <- as.matrix(statMixRHLP$tau_ijgk[, , g]) # [(nxm) x R]

        if (variance_type == "homoskedastic") {
          s <- 0
        } else {
          sigma_gk <- zeros(R, 1)
        }

        beta_gk <- matrix(NA, p + 1, R)

        for (k in 1:R) {

          segments_weights <- tauijk[, k, drop = F]
          phigk <- (sqrt(cluster_weights * segments_weights) %*% ones(1, p + 1)) * phi$XBeta #[(n*m)*(p+1)]
          Xgk <- sqrt(cluster_weights * segments_weights) * fData$vecY

          # Maximization w.r.t beta_gk: Weighted least squares
          beta_gk[, k] <- solve(t(phigk) %*% phigk + .Machine$double.eps * diag(p + 1)) %*% t(phigk) %*% Xgk # Maximization w.r.t betagk

          #    the same as
          #                 W_gk = diag(cluster_weights.*segment_weights);
          #                 beta_gk(:,k) = inv(phiBeta'*W_gk*phiBeta)*phiBeta'*W_gk*X;
          #   Maximization w.r.t au sigma_gk :
          if (variance_type == "homoskedastic") {
            sk <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2)
            s <- s + sk
            sigma_gk <- s / sum(colSums((cluster_weights %*% ones(1, R)) * tauijk))
          } else {
            sigma_gk[k] <- colSums((Xgk - phigk %*% beta_gk[, k]) ^ 2) / (colSums(cluster_weights * segments_weights))
          }
        }

        beta[, , g] <<- beta_gk
        if (variance_type == "homoskedastic") {
          sigma2[g] <<- sigma_gk
        } else {
          sigma2[, g] <<- sigma_gk
        }

        # Maximization w.r.t W
        # Setting of W[,,g]
        Wg_init <- matrix(W[, , g], nrow = q + 1)

        res_irls <- IRLS(phi$Xw, tauijk, cluster_weights, Wg_init, verbose_IRLS)

        W[, , g] <<- res_irls$W
        piik <- res_irls$piik

      }
    }
  )
)
