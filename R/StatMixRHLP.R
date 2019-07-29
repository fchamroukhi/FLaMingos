#' A Reference Class which contains statistics of a MixRHLP model.
#'
#' StatMixRHLP contains all the statistics associated to a
#' [MixRHLP][ParamMixRHLP] model.
#'
#' @field pi_jgk Array of size \eqn{(nm, R, K)} representing the logistic
#'   proportion for cluster g.
#' @field h_ig Matrix of size \eqn{(n, K)} giving the posterior probabilities
#'   (fuzzy segmentation matrix) that the curve \eqn{Y_{i}} originates from the
#'   \eqn{g}-th RHLP model.
#' @field c_ig Hard segmentation logical matrix of dimension \eqn{(n, K)}
#'   obtained by the Maximum a posteriori (MAP) rule: \eqn{c\_ig = 1 \
#'   \textrm{if} \ c\_ig = \textrm{arg} \ \textrm{max}_{s} \ h\_is;\ 0 \
#'   \textrm{otherwise}}{c_ig = 1 if c_ig = arg max_s h_is; 0 otherwise}, \eqn{g
#'   = 1,\dots,K}.
#' @field klas Column matrix of the labels issued from `c_ig`. Its elements are
#'   \eqn{klas(i) = g}, \eqn{i = 1,\dots,n}.
#' @field tau_ijgk Array of size \eqn{(nm, R, K)} giving the posterior
#'   probabilities that the observation \eqn{Y_{ij}} originates from the
#'   \eqn{R}-th regime of the \eqn{g}-th RHLP model.
#' @field polynomials Array of size \eqn{(m, R, K)} giving the values of the
#'   estimated polynomial regression components.
#' @field weighted_polynomials Array of size \eqn{(m, R, K)} giving the values
#'   of the estimated polynomial regression components weighted by the prior
#'   probabilities `pi_jgk`.
#' @field Ex Matrix of size \emph{(m, K)}. `Ex` is the curve expectation
#'   (estimated signal): sum of the polynomial components weighted by the
#'   logistic probabilities `pi_jgk`.
#' @field loglik Numeric. Observed-data log-likelihood of the MixRHLP model.
#' @field com_loglik Numeric. Complete-data log-likelihood of the MixRHLP model.
#' @field stored_loglik Numeric vector. Stored values of the log-likelihood at
#'   each EM iteration.
#' @field stored_com_loglik Numeric vector. Stored values of the Complete
#'   log-likelihood at each EM iteration.
#' @field BIC Numeric. Value of BIC (Bayesian Information Criterion).
#' @field ICL Numeric. Value of ICL (Integrated Completed Likelihood).
#' @field AIC Numeric. Value of AIC (Akaike Information Criterion).
#' @field log_fg_xij Matrix of size \eqn{(n, K)} giving the values of the
#'   probability density function \eqn{f(y_{i} | h_ig = 1, \boldsymbol{x},
#'   \boldsymbol{\Psi})}{f(y_{i} | h_ig = 1, x, \Psi)}, \eqn{i = 1,\dots,n}.
#' @field log_alphag_fg_xij Matrix of size \eqn{(n, K)} giving the values of the
#'   logarithm of the joint probability density function \eqn{f(y_{i}, \ h_{i} =
#'   g | \boldsymbol{x}, \boldsymbol{\Psi})}{f(y_{i}, h_{i} = g | x, \Psi)},
#'   \eqn{i = 1,\dots,n}.
#' @field log_tau_ijgk Array of size \eqn{(nm, R, K)} giving the logarithm of
#'   `tau_ijgk`.
#' @seealso [ParamMixRHLP]
#' @export
StatMixRHLP <- setRefClass(
  "StatMixRHLP",
  fields = list(
    pi_jgk = "array", # pi_jgk :logistic proportions for cluster g
    h_ig = "matrix", # h_ig = prob(curve|cluster_g) : post prob (fuzzy segmentation matrix of dim [nxK])
    c_ig = "matrix", # c_ig : Hard partition obtained by the MAP rule:  c_{ig} = 1
    # if and only c_i = arg max_g h_ig (g=1,...,K)
    klas = "matrix", # klas : column vector of cluster labels
    Ex = "matrix",
    # Ex: curve expectation: sum of the polynomial components beta_gk weighted by
    # the logitic probabilities pi_jgk: Ex(j) = sum_{k=1}^R pi_jgk beta_gk rj, j=1,...,m. Ex
    # is a column vector of dimension m for each g.
    loglik = "numeric", # the loglikelihood of the EM or CEM algorithm
    com_loglik = "numeric", # the complete loglikelihood of the EM (computed at the convergence) or CEM algorithm
    stored_loglik = "numeric", # vector of stored valued of the comp-log-lik at each EM teration
    stored_com_loglik = "numeric",
    tau_ijgk = "array",
    # tau_ijgk prob(y_{ij}|kth_segment,cluster_g), fuzzy
    # segmentation for the cluster g. matrix of dimension
    # [nmxR] for each g  (g=1,...,K).
    log_tau_ijgk = "array",
    BIC = "numeric", # BIC value = loglik - nu*log(nm)/2.
    ICL = "numeric", # ICL value = comp-loglik_star - nu*log(nm)/2.
    AIC = "numeric", # AIC value = loglik - nu.
    log_fg_xij = "matrix",
    log_alphag_fg_xij = "matrix",
    polynomials = "array",
    weighted_polynomials = "array"
  ),
  methods = list(
    initialize = function(paramMixRHLP = ParamMixRHLP()) {

      pi_jgk <<- array(0, dim = c(paramMixRHLP$fData$m * paramMixRHLP$fData$n, paramMixRHLP$R, paramMixRHLP$K))
      h_ig <<- matrix(NA, paramMixRHLP$fData$n, paramMixRHLP$K)
      c_ig <<- matrix(NA, paramMixRHLP$fData$n, paramMixRHLP$K)
      klas <<- matrix(NA, paramMixRHLP$fData$n, 1)
      Ex <<- matrix(NA, paramMixRHLP$fData$m, paramMixRHLP$K)
      loglik <<- -Inf
      com_loglik <<- -Inf
      stored_loglik <<- numeric()
      stored_com_loglik <<- numeric()
      BIC <<- -Inf
      ICL <<- -Inf
      AIC <<- -Inf
      log_fg_xij <<- matrix(0, paramMixRHLP$fData$n, paramMixRHLP$K)
      log_alphag_fg_xij <<- matrix(0, paramMixRHLP$fData$n, paramMixRHLP$K)
      polynomials <<- array(NA, dim = c(paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
      weighted_polynomials <<- array(NA, dim = c(paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
      tau_ijgk <<- array(0, dim = c(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
      log_tau_ijgk <<- array(0, dim = c(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$R, paramMixRHLP$K))
    },

    MAP = function() {
      "MAP calculates values of the fields \\code{c_ig} and \\code{klas}
      by applying the Maximum A Posteriori Bayes allocation rule.

      \\eqn{c\\_ig = 1 \\ \\textrm{if} \\ c\\_ig = \\textrm{arg} \\
      \\textrm{max}_{s} \\ h\\_is;\\ 0 \\ \\textrm{otherwise}
      }{c_ig = 1 if c_ig = arg max_s h_is; 0 otherwise}, \\eqn{g = 1,\\dots,K}."

      N <- nrow(h_ig)
      R <- ncol(h_ig)
      ikmax <- max.col(h_ig)
      ikmax <- matrix(ikmax, ncol = 1)
      c_ig <<- ikmax %*% ones(1, R) == ones(N, 1) %*% (1:R)
      klas <<- ones(N, 1)
      for (k in 1:R) {
        klas[c_ig[, k] == 1] <<- k
      }
    },

    computeStats = function(paramMixRHLP) {
      "Method used in the EM algorithm to compute statistics based on
      parameters provided by the object \\code{paramMixRHLP} of class
      \\link{ParamMixRHLP}."

      for (g in 1:paramMixRHLP$K) {

        polynomials[, , g] <<- paramMixRHLP$phi$XBeta[1:paramMixRHLP$fData$m, ] %*% matrix(paramMixRHLP$beta[, , g], nrow = paramMixRHLP$p + 1)

        weighted_polynomials[, , g] <<- pi_jgk[1:paramMixRHLP$fData$m, , g] * polynomials[, , g]
        Ex[, g] <<- rowSums(as.matrix(weighted_polynomials[, , g]))

      }

      Ex <<- matrix(Ex, nrow = paramMixRHLP$fData$m)

      BIC <<- loglik - (paramMixRHLP$nu * log(paramMixRHLP$fData$n) / 2)
      AIC <<- loglik - paramMixRHLP$nu

      cig_log_alphag_fg_xij <- (c_ig) * (log_alphag_fg_xij)

      com_loglik <<- sum(rowSums(cig_log_alphag_fg_xij))

      ICL <<- com_loglik - paramMixRHLP$nu * log(paramMixRHLP$fData$n) / 2
    },

    CStep = function(reg_irls) {
      "Method used in the CEM algorithm to update statistics."

      h_ig <<- exp(lognormalize(log_alphag_fg_xij))

      MAP() # Setting klas and c_ig

      # Compute the optimized criterion
      cig_log_alphag_fg_xij <- (c_ig) * log_alphag_fg_xij
      com_loglik <<- sum(cig_log_alphag_fg_xij) +  reg_irls
    },

    EStep = function(paramMixRHLP) {
      "Method used in the EM algorithm to update statistics based on parameters
      provided by the object \\code{paramMixRHLP} of class \\link{ParamMixRHLP}
      (prior and posterior probabilities)."

      for (g in 1:paramMixRHLP$K) {

        alpha_g <- paramMixRHLP$alpha[g]
        beta_g <- matrix(paramMixRHLP$beta[, , g], nrow = paramMixRHLP$p + 1)
        Wg <- matrix(paramMixRHLP$W[, , g], nrow = paramMixRHLP$q + 1)
        piik <- multinomialLogit(Wg, paramMixRHLP$phi$Xw, ones(nrow(paramMixRHLP$phi$Xw), ncol(Wg) + 1), ones(nrow(paramMixRHLP$phi$Xw), 1))$piik
        pi_jgk[, , g] <<- as.matrix(repmat(piik[1:paramMixRHLP$fData$m,], paramMixRHLP$fData$n, 1))

        log_pijgk_fgk_xij <- zeros(paramMixRHLP$fData$n * paramMixRHLP$fData$m, paramMixRHLP$R)

        for (k in 1:paramMixRHLP$R) {

          beta_gk <- as.matrix(beta_g[, k])
          if (paramMixRHLP$variance_type == "homoskedastic") {
            sgk <- paramMixRHLP$sigma2[g]
          } else {
            sgk <- paramMixRHLP$sigma2[k, g]
          }
          z <- ((paramMixRHLP$fData$vecY - paramMixRHLP$phi$XBeta %*% beta_gk) ^ 2) / sgk
          log_pijgk_fgk_xij[, k] <- log(pi_jgk[, k, g]) - 0.5 * (log(2 * pi) + log(sgk)) - 0.5 * z # pdf cond c_i = g et z_i = k de xij
        }

        log_pijgk_fgk_xij <- pmin(log_pijgk_fgk_xij, log(.Machine$double.xmax))
        log_pijgk_fgk_xij <- pmax(log_pijgk_fgk_xij, log(.Machine$double.xmin))

        pijgk_fgk_xij <- exp(log_pijgk_fgk_xij)
        sumk_pijgk_fgk_xij <- rowSums(pijgk_fgk_xij) # sum over k
        log_sumk_pijgk_fgk_xij <- log(sumk_pijgk_fgk_xij) # [n*m, 1]

        log_tau_ijgk[, , g] <<- log_pijgk_fgk_xij - log_sumk_pijgk_fgk_xij %*% ones(1, paramMixRHLP$R)
        tau_ijgk[, , g] <<- exp(lognormalize(log_tau_ijgk[, , g]))

        log_fg_xij[, g] <<- rowSums(t(matrix(log_sumk_pijgk_fgk_xij, paramMixRHLP$fData$m, paramMixRHLP$fData$n))) # [n x 1]:  sum over j=1,...,m: fg_xij = prod_j sum_k pi_{jgk} N(x_{ij},mu_{gk},s_{gk))
        log_alphag_fg_xij[, g] <<- log(alpha_g) + log_fg_xij[, g] # [nxg]
      }
      log_alphag_fg_xij <<- pmin(log_alphag_fg_xij, log(.Machine$double.xmax))
      log_alphag_fg_xij <<- pmax(log_alphag_fg_xij, log(.Machine$double.xmin))

      h_ig <<- exp(lognormalize(log_alphag_fg_xij))

      loglik <<- sum(log(rowSums(exp(log_alphag_fg_xij))))
    }

  )
)
