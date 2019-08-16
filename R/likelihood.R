#' @importFrom edgeR dispCoxReid

# This is adopted from a DESeq function
profileLogLikelihood2 <- function( disp, mm, y, muhat, type = "CR",
                                   useCpp = F)
{
  #tic = proc.time()

  # calculate the log likelihood:
  if(length(disp) != length(y)){
    disp <- rep(disp, length(y))
  }

  ## !!! wrong Cpp code when 1 / disp is very large, e.g., 1e30
  ## This might be a problem in optimization when disp is very small
  ## So R function dnbinom gives the accurate result
  if (!useCpp) {
    ll = sum(dnbinom(y, mu = muhat, size = 1 / disp, log = T))
  } else {
    ll = likelihoodNegativeBinomial(y, muhat, 1/disp)
  }
  #print(system.time({
  # ll <- sum( sapply( seq(along=y), function(i)
  # dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) ) )
  #dnbinomLog(y[i], muhat[i], 1/disp[i])))
  #}))

  # toc = proc.time()
  # print(toc - tic)

  # transform the residuals, i.e., y - muhat, to the linear
  # predictor scale by multiplying them with the derivative
  # of the link function, i.e., by 1/muhat, and add this to the
  # linear predictors, log(muhat), to get the predictors that
  # are used in the IWLS regression
  # z <- log(muhat) + ( y - muhat ) / muhat

  if (type == "CR") {
    #tic = proc.time()

    # the variance function of the NB is as follows
    #print(system.time({
    v0 <- muhat + disp * muhat^2

    # transform the variance vector to linear predictor scale by
    # multiplying with the squared derivative of the link to
    # get the (reciprocal) weights for the IWLS
    sqrtw <- sqrt(abs(muhat) / (1 + muhat * disp))
    #}))


    # All we need from the IRLS run is the QR decomposition of
    # its matrix
    #print(system.time({
    qrres <- qr( mm* c(sqrtw) )
    #}))


    # from it, we extract we leverages and calculate the Cox-Reid
    # term:
    #print(system.time({
    cr <- sum( log( abs( diag( qrres$qr )[ seq_len(qrres$rank) ] ) ) )
    #}))

    # return the profile log likelihood:
    ll = ll - cr
  }

  # toc = proc.time()
  # print(toc - tic)
  return(ll)
}

# This is adopted from a DESeq funciton
estimateAndFitDispersionsWithCoxReid <- function( counts, modelFormula, modelFrame,
                                                  sizeFactors, fitType = c( "parametric", "local" ),
                                                  locfit_extra_args=list(), lp_extra_args=list(), initialGuess=.1,
                                                  dispMethod = "poisson-CR")
{
  # browser()
  if( as.character(modelFormula[2]) == "count" )
    modelFormula <- modelFormula[-2]   # the '[-2]' removes the lhs, i.e., the response
  # browser()
  if (is.data.frame(modelFrame) && length(unique(modelFrame[, 1])) == 1) {
    # check the first column
    mm = rep(1, dim(counts)[2])
  } else {
    mm <- model.matrix( modelFormula, modelFrame )
  }
  disps <- apply( counts, 1, function( y ) {
    #browser()
    if (dispMethod == "pooled-CR" || dispMethod == "poisson-CR" ||
        dispMethod == "pooled-ML" || dispMethod == "poisson-ML") {
      type = "CR"
      if (dispMethod == "pooled-ML" || dispMethod == "poisson-ML") {
        type = "ML"
      }
      fit <- try({
        #tic = proc.time()

        if (dispMethod == "pooled-CR" || dispMethod == "pooled-ML") {
          fitted = suppressWarnings(glm.fit( mm, y,
                              family=MASS::negative.binomial( initialGuess ),
                              offset=log(sizeFactors) ))
        } else if (dispMethod == "poisson-CR" || dispMethod == "poisson-ML") {
          fitted = glm.fit( mm, y, family=poisson(),
                            #control = glm.control(epsilon = 1e-4, maxit = 25, trace = FALSE),
                            offset=log(sizeFactors))
        }

        #toc = proc.time()
        #print(toc - tic)
        #fitted
      },
      silent=TRUE)

      if( inherits( fit, "try-error" ) )
        NA
      else {
        if( df.residual(fit) == 0 ) {
          stop( "No residual degrees of freedom. Most likely the design is lacking sufficient replication." )
        }

        #browser()

        # # linear search
        # muhat = fitted.values(fit)
        # tic = proc.time()
        # optimizeResult = optimize(
        #   function(logalpha)
        #     profileLogLikelihood( exp(logalpha), mm, y, muhat ),
        #   #log( c( 1e-11, 1e5 ) ),
        #   c(-16.1206, -16.1208),
        #   maximum=TRUE
        # )
        # value = exp(optimizeResult$maximum)
        # toc = proc.time()
        # print("optimize")
        # print(toc - tic)
        # print(optimizeResult)

        #tic = proc.time()
        muhat = fitted.values(fit)
        optimizeResult = optimize(
          function(logalpha)
            profileLogLikelihood2( exp(logalpha), mm, y, muhat, type),
          log( c( 1e-11, 1e5 ) ),
          #c(0.44372, 0.44373),
          maximum=TRUE
        )
        value = exp(optimizeResult$maximum)
        # print(value)
        #toc = proc.time()
        #print("optimize")
        #print(toc - tic)
        # print(optimizeResult)


        # # general purpose opt
        # tic = proc.time()
        # # optimResult = optim(0, function(logalpha)
        # #   profileLogLikelihood( exp(logalpha), mm, y, fitted.values(fit) ),
        # #   method = "Brent",
        # #   lower = log(1e-11), upper = log(1e5),
        # #   control = list(fnscale = -1))
        # # value = exp(optimResult$par)
        #
        # optimResult = optimize(
        #   function(logalpha)
        #     profileLogLikelihood( exp(logalpha), mm, y, fitted.values(fit) ),
        #   log( c( 0.01, 64 ) ),
        #   maximum=TRUE
        # )
        # value = exp(optimizeResult$maximum)
        # toc = proc.time()
        # print("optim")
        # print(toc - tic)
        # print(optimResult)

        return(value)
      }
    } else {
      # use edgeR's function
      if (length(unique(modelFrame[, 1])) == 1) {
        design = NULL
      } else {
        design = mm
      }
      disps = dispCoxReid(matrix(y, nrow=1), design = design,
                          interval = c( 1e-11, 1e5 ),
                          offset=log(sizeFactors),
                          min.row.sum=0)
      return(disps)
    }
  }
  )

  means <- colMeans( t(counts) / sizeFactors )
  xim <- mean( 1/sizeFactors )

  ans = NA
  if (dispMethod == "pooled-CR" || dispMethod == "full-CR" ||
      dispMethod == "poisson-CR" || dispMethod == "pooled-ML" ||
      dispMethod == "poisson-ML") {

  } else if( fitType == "local" ) {

    fit <- do.call( "locfit", c(
      list(
        disps[means>0] ~ do.call( "lp", c( list( log(means[means>0]) ), lp_extra_args ) ),
        family = "gamma" ),
      locfit_extra_args ) )

    rm( means )

    ans <- function( q )
      pmax( ( safepredict( fit, log(q) ) - xim * q ) / q^2, 1e-8 )

    # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
    # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.

  } else if( fitType == "parametric" )

    ans <- parametricDispersionFit( means, disps )

  else
    stop( "Unkknown fitType." )

  attr( ans, "fitType" ) <- fitType
  list( disps=disps, dispFunc=ans )

}

NBDispersionAndLikelihood = function(count, modelFormula, modelFrame,
                                     countPerCell, dispMethod,
                                     groups = NULL) {
  # assume no covariates when using glm.nb

  # browser()
  dispersion = NA
  logLik = NA
  seTheta = NA
  if (dispMethod == "full-CR" || dispMethod == "pooled-CR" ||
      dispMethod == "poisson-CR" || dispMethod == "pooled-ML" ||
      dispMethod == "poisson-ML") {
    dispersionPrevious = 0.1
    logLikPrevious = NULL
    iter = 0
    iterMax = 50
    eps = 1e-2
    stepsize = 1e-2
    converged = F
    while (iter <= iterMax) {
      #browser()
      dispersionNew = estimateAndFitDispersionsWithCoxReid(
        matrix(count, nrow = 1),
        modelFormula, modelFrame,
        countPerCell, "local", list(), list(),
        1 / dispersionPrevious,
        dispMethod)$disps
      if (iter == 0) {
        dispersion = dispersionNew
      } else {
        dispersion = dispersionPrevious + (dispersionNew - dispersionPrevious) *
          stepsize
      }
      # # <debug>
      # dispersion = 1
      # # </debug>
      if( as.character(modelFormula[2]) == "count" ) {
        modelFormula <- modelFormula[-2]
      }
      mm <- model.matrix( modelFormula, modelFrame )
      model = suppressWarnings(glm.fit( mm, count,
                       family=MASS::negative.binomial(1 / dispersion),
                       offset=log(countPerCell),
                       control = glm.control(maxit = 30)))
      # # <debug>
      # model = glm.fit( c(1, 1, 2, 2), c(1, 2, 3, 4),
      #                  family=poisson(),
      #                  offset=log(countPerCell[1:4]) )
      # # </debug>
      logLik = -0.5 * (model$aic - 2 * (length(count) - model$df.residual))

      if (!is.null(logLikPrevious) && (logLik - logLikPrevious < eps)) {
        converged = modelPrevious$converged
        break
      }

      logLikPrevious = logLik
      dispersionPrevious = dispersion
      modelPrevious = model
      iter = iter + 1
    }
  } else if (dispMethod == "ML") {
    try({
      model = suppressWarnings(glm.nb(formula(paste("count ~ groups",
                                   "+ offset(log(countPerCell))"))))
      aicNB = model$aic
      logLik = -0.5 * (model$aic - 2 * 3)
      dispersion = 1 / model$theta
      seTheta = model$SE.theta
    })
  } else {
    stop("wrong dispersion estimation method")
  }
  # glmFit(matrix(count, nrow = 1), mm, dispersion = dispersion,
  #                 prior.count = 0)
  return(list(dispersion = dispersion, logLik = logLik, model = model,
              seTheta = seTheta, converged = converged))
}

