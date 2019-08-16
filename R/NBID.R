#' @import pscl
#' @import MASS
#' @import nloptr

dispersionLowerBound = 1e-11
dispersionUpperBound = 1e5

NBGeneral <- function(r) {
  # modified from deseq and MASS

  # TODO: need to make sure it is correct with multiple dispersions
  fam <- list(
    family = sprintf( "NBGeneral, r = 1 / dispersion and is a vector"),
    link = "log",
    linkfun = function(mu) log( mu),
    linkinv = function (eta) pmax(exp(eta), .Machine$double.eps),
    mu.eta = function (eta) pmax(exp(eta), .Machine$double.eps),
    variance = function(mu) mu + mu^2 / r,
    dev.resids = function(y, mu, wt) {
      2 * wt * ( ifelse( y > 0, y * log( y / mu ), 0 ) +
                   (r + y) * log( (r+mu)/(r+y) ) )
      },
    aic = function(y, n, mu, wt, dev) {   # it is actually -2 * logLik
      term <- (y + r) * log(mu + r) - y * log(mu) +
        lgamma(y + 1) - r * log(r) + lgamma(r) - lgamma(r+y)
      2 * sum(term * wt)
    },
    initialize = expression( {
      n <- rep.int(1, nobs)
      mustart <- y + 0.1
    } ),
    valid.mu <- function(mu) all( mu > 0 ),
    valid.eta <- function(eta) TRUE,
    simulate <- NA
  )

  class(fam) <- "family"
  return(fam)
}

likelihoodDifferentDispersionNegative = function(parameters, count,
                                       groups, modelMatrix, countPerCell,
                                       dispersionsGroupName,
                                       type) {
  nDispersion = length(dispersionsGroupName)
  parameters[1:nDispersion] = exp(parameters[1:nDispersion])
  value = -1 * likelihoodDifferentDispersion(parameters, count,
                                        groups, modelMatrix, countPerCell,
                                        dispersionsGroupName,
                                        type)
  return(value)
}

likelihoodDifferentLogDispersion = function(parameters, otherData) {
  count = otherData$count
  groups = otherData$groups
  modelMatrix = otherData$modelMatrix
  countPerCell = otherData$countPerCell
  dispersionsGroupName = otherData$dispersionsGroupName
  type = "ML"
  nDispersion = length(dispersionsGroupName)

  # fill in the same dispersion for all groups
  # browser()
  if (nDispersion > 1 && length(parameters) - dim(modelMatrix)[2] == 1) {
    parameters = c(rep(parameters[1], nDispersion - 1), parameters)
  } else if (nDispersion == 0 && length(parameters) == dim(modelMatrix)[2]) {
    nGroup = length(unique(groups))
    parameters = c(rep(otherData$dispersion, nGroup),
                   parameters)
    dispersionsGroupName = unique(groups)
    nDispersion = nGroup
  }

  if (nDispersion > 0) {
    parameters[1:nDispersion] = exp(parameters[1:nDispersion])
  }
  value = likelihoodDifferentDispersion(parameters, count,
                                             groups, modelMatrix, countPerCell,
                                             dispersionsGroupName,
                                             type)
  return(value)
}

likelihoodDifferentLogDispersionNegative = function(parameters, otherData) {
  negativeValue = - likelihoodDifferentLogDispersion(parameters, otherData)
  return(negativeValue)
}

calculateDispersionDerivative = function(c, y, mu) {
  firstPart = 1 / c * sum(log(1 + c * mu))
  secondPart = sum(mu * (1 + c * y) / (1 + c * mu))
  thirdPart = 0
  for (i in 1:length(y)) {
    if (y[i] >= 2) {
      s = c * (1:(y[i] - 1))
      thirdPart = thirdPart + sum(s / (1 + s))
    }
  }

  derivative = firstPart - secondPart + thirdPart
  if (is.na(derivative)) {
    print("dump")
    print(c)
    print(firstPart)
    print(secondPart)
    print(thirdPart)
    print("dump end")
  }
  return(derivative)
}

gradientDifferentLogDispersion = function(parameters, otherData) {
  # calculate the gradient value

  count = otherData$count
  groups = otherData$groups
  modelMatrix = otherData$modelMatrix
  countPerCell = otherData$countPerCell
  dispersionsGroupName = otherData$dispersionsGroupName


  nDispersion = length(dispersionsGroupName)
  dispersionVector = numeric(length(groups))
  dispersionVector[] = dispersionLowerBound
  if (nDispersion > 0) {
    parameters[1:nDispersion] = exp(parameters[1:nDispersion])
    for (i in 1:nDispersion) {
      dispersionVector[groups == dispersionsGroupName[i]] = parameters[i]
    }
    beta = parameters[-(1:nDispersion)]
  } else {
    beta = parameters
  }

  #browser()
  mu = c(countPerCell * exp(modelMatrix %*% beta))
  v = mu * (1 + mu * dispersionVector)
  D = modelMatrix * mu

  # gradient of beta
  gradBeta = t(D) %*% ((1 / v) * (count - mu))

  # browser()
  # gradient with respect to phi = log(c), c is the dispersion
  gradDispersion = numeric(nDispersion)
  if (nDispersion > 0) {
    for (i in 1:nDispersion) {
      index = which(groups == dispersionsGroupName[i])
      c = parameters[i]
      gradDispersion[i] = calculateDispersionDerivative(c, count[index],
                                                        mu[index])
    }
  }

  gradient = c(gradDispersion, gradBeta)

  return(gradient)
}

gradientDifferentLogDispersionNegative = function(parameters, otherData) {
  negativeValues = -gradientDifferentLogDispersion(parameters, otherData)
  return(negativeValues)
}

calculateDispersionSecondDerivative = function(c, y, mu) {
  # input
  # c: the dispersion, a single value
  # y: the count
  # mu: the mean

  secondDerivative = - 1 / c * sum(log(1 + c * mu)) + sum(mu / (1 + c * mu)) +
                     - c * sum(mu * (y - mu) / (1 + c * mu)^2)
  for (i in 1:length(y)) {
    if (y[i] >= 2) {
      s = 1:(y[i] - 1)
      secondDerivative = secondDerivative + sum(s / (1 + s)^2)
    }
  }

  return(secondDerivative)
}

calculateDispersionBetaSecondDerivative = function(c, y, mu, modelMatrix) {
  # input
  # c: the dispersion, a single value

  value = - c * (y - mu) * mu / (1 + c * mu)^2
  secondDerivative = t(modelMatrix) %*% value
  return(secondDerivative)
}

calculateBetaSecondDerivative = function(cVector, y, mu, X) {
  # input
  # cVector: dispersion, a vector with the same length of y and mu

  w = mu * (1 + cVector * y) / (1 + cVector * mu)^2
  # H = - Xt %*% w %*% X
  hessianMatrix = -t(X) %*% (X * w)
  return(hessianMatrix)
}

hessianDifferentLogDispersion = function(parameters, otherData) {
  # the second derivatives

  #browser()

  # set up the data
  count = otherData$count
  groups = otherData$groups
  modelMatrix = otherData$modelMatrix
  countPerCell = otherData$countPerCell
  dispersionsGroupName = otherData$dispersionsGroupName

  nDispersion = length(dispersionsGroupName)
  dispersionVector = numeric(length(groups))
  dispersionVector[] = dispersionLowerBound
  if (nDispersion > 0) {
    parameters[1:nDispersion] = exp(parameters[1:nDispersion])
    for (i in 1:nDispersion) {
      dispersionVector[groups == dispersionsGroupName[i]] = parameters[i]
    }
    beta = parameters[-(1:nDispersion)]
  } else {
    beta = parameters
  }

  mu = c(countPerCell * exp(modelMatrix %*% beta))

  # calculate the hessian matrix
  nParameters = length(parameters)
  hessianMatrix = matrix(NA, length(parameters), nParameters)
  if (nDispersion > 0) {
    hessianDispersion = numeric(nDispersion)
    hessianDispersionBeta = matrix(NA, nDispersion, nParameters - nDispersion)
    for (i in 1:nDispersion) {
      # fill in the dispersion block
      index = which(groups == dispersionsGroupName[i])
      c = parameters[i]
      hessianDispersion[i] = calculateDispersionSecondDerivative(
                          c, count[index],
                          mu[index]
      )

      # fill in the dispersion-beta block
      hessianDispersionBeta[i, ] = calculateDispersionBetaSecondDerivative(
                                     c, count[index],
                                     mu[index], modelMatrix[index, , drop = F]
      )
    }
    hessianMatrix[1:nDispersion, 1:nDispersion] = diag(hessianDispersion)
    hessianMatrix[1:nDispersion, (nDispersion + 1):nParameters] =
      hessianDispersionBeta
    hessianMatrix[(nDispersion + 1):nParameters, 1:nDispersion] =
      t(hessianDispersionBeta)
  }

  # browser()
  # fill in the beta block
  hessianBeta = calculateBetaSecondDerivative(dispersionVector, count, mu,
                                              modelMatrix)
  hessianMatrix[(nDispersion + 1):nParameters, (nDispersion + 1):nParameters] =
    hessianBeta

  return(hessianMatrix)
}

likelihoodDifferentDispersion = function(parameters, count,
                                         groups, modelMatrix, countPerCell,
                                         dispersionsGroupName,
                                         type) {
  # calculate the likelihood of the data based on parameters
  # the first two are dispersions, the rest are beta for the linear predictors

  nDispersion = length(dispersionsGroupName)
  dispersionVector = numeric(length(groups))
  dispersionVector[] = dispersionLowerBound
  if (nDispersion > 0) {
    for (i in 1:nDispersion) {
      dispersionVector[groups == dispersionsGroupName[i]] = parameters[i]
    }
    beta = parameters[-(1:nDispersion)]
  } else {
    beta = parameters
  }

  #browser()
  mu = countPerCell * exp(modelMatrix %*% beta)
  logLikelihood = profileLogLikelihood2(dispersionVector, modelMatrix,
                                        count, mu,
                                        type = type)
  return(logLikelihood)
}

likelihoodWithKnownMean = function(dispersion, mu,
                                   count, groups, modelMatrix,
                                   type) {
  # input
  # dispersion: a vector of two elements corresponds to two group factors
  # groups: a factor

  dispersionVector = numeric(length(groups))
  dispersionVector[groups == 1] = dispersion[1]
  dispersionVector[groups == 2] = dispersion[2]

  logLikelihood = profileLogLikelihood2(dispersionVector, modelMatrix,
                                        count, mu,
                                        type = type)
  return(logLikelihood)
}

MLENBDirect = function(count, groups, covariates, countPerCell,
                               dispersions, dispersionsGroupName,
                             type = "ML") {
  dispersionsOptim = NA
  logLikOptim = NA
  coef = NA

  groups = factor(groups)
  if (!is.null(covariates)) {
    modelMatrix = model.matrix(~covariates)
  } else {
    modelMatrix = matrix(1, length(count), 1)
  }

  dispersionVector = numeric(length(groups))
  dispersionVector[] = dispersionLowerBound
  nDispersion = length(dispersions)
  if (nDispersion > 0) {
    for (i in 1:nDispersion) {
      dispersionVector[groups == dispersionsGroupName[i]] = dispersions[i]
    }
  }

  modelInit = suppressWarnings(glm.fit(modelMatrix, count,
                     #family = NBGeneral(dispersionVector),
                     family = MASS::negative.binomial( 1 / dispersionVector),
                     offset = log(countPerCell)))
  parameters = c(dispersions, coefficients(modelInit))

  # browser()
  # result = optim(parameters, likelihoodDifferentDispersion, gr = NULL,
  #                count = count, groups = groups,
  #                modelMatrix = modelMatrix,
  #                countPerCell = countPerCell,
  #                type = type,
  #                #method = "Nelder-Mead",
  #                #lower = 1e-7, upper = 1e5,
  #                control = list(fnscale = -1))
  # dispersions = result$par[1:nDispersion]
  # logLikOptim = result$value
  # coef = result$par[-(1:nDispersion)]

  # use nloptr
  # # do a globle search
  # opts = list("algorithm"="NLOPT_GN_DIRECT_L",
  #             "xtol_rel"=1.0e-8,
  #             "maxeval" = 1e3)
  # result = nloptr(c(log(parameters[1]), parameters[2:length(parameters)]),
  #                 eval_f = likelihoodDifferentDispersionNegative,
  #                 lb = c(-11, rep(-1e5, dim(modelMatrix)[2])),
  #                 ub = c(5, rep(1e5, dim(modelMatrix)[2])),
  #                 opts = opts,
  #                 count = count, groups = groups,
  #                 modelMatrix = modelMatrix,
  #                 countPerCell = countPerCell,
  #                 type = type)
  # parameters = result$solution

  if (nDispersion > 0) {
    # do a local search
    opts = list("algorithm"="NLOPT_LN_NEWUOA", #"NLOPT_LN_NELDERMEAD",
                 "xtol_rel"=1.0e-8,
                "maxeval" = 1e4)
    result = nloptr::nloptr(c(log(parameters[1:nDispersion]),
                      parameters[(nDispersion + 1):length(parameters)]),
                    eval_f = likelihoodDifferentDispersionNegative,
                    lb = c(rep(log(dispersionLowerBound), nDispersion),
                           rep(-Inf, length(parameters) - nDispersion)),
                    ub = c(rep(log(dispersionUpperBound), nDispersion),
                           rep(Inf, length(parameters) - nDispersion)),
                    opts = opts,
                    count = count, groups = groups,
                    modelMatrix = modelMatrix,
                    countPerCell = countPerCell,
                    dispersionsGroupName = dispersionsGroupName,
                    type = type)
    #browser()
    if (result$status >= 0 && result$status <= 4) {
      dispersionsOptim = exp(result$solution[1:nDispersion])
      logLikOptim = -result$objective
      coef = result$solution[-(1:nDispersion)]
    } else {
      print("nloptr failed!")
      browser()
    }
  } else {
    # browser()
    dispersionsOptim = NULL
    logLikOptim = -0.5 * (modelInit$aic - 2 * dim(modelMatrix)[2])
    coef = coef(modelInit)
  }

  # # <debug>
  # dispersionVector = numeric(length(groups))
  # for (i in 1:nDispersion) {
  #   dispersionVector[groups == i] = dispersions[i]
  # }
  # model = suppressWarnings(glm.fit(modelMatrix, count,
  #                    #family = NBGeneral(dispersionVector),
  #                    family = MASS::negative.binomial( 1 / dispersionVector),
  #                    offset = log(countPerCell)))
  # mu = model$fitted.values
  # logLik = -0.5 * (model$aic - 2 * 1) # assume no covariates
  # # </debug>

  return(list(dispersions = dispersionsOptim, logLik = logLikOptim,
              coef = coef))
}

MLEUnderH0EqualMean = function(count, groups, covariates, countPerCell,
                               dispersions) {
  # assume two groups

  # input
  # dispersions: the inital values of dispersions

  #browser()
  groups = factor(groups)
  if (!is.null(covariates)) {
    modelMatrix = model.matrix(~covariates)
  } else {
    modelMatrix = matrix(1, length(count), 1)
  }

  # # use Poisson distribution to estimate the mean
  # modelInit = glm.fit(modelMatrix, count, family = poisson(),
  #                     offset = log(countPerCell))
  # mu = modelInit$fitted.values

  logLikPrevious = NULL
  dispersionsPrevious = NULL
  eps = 1e-4
  epsDispersion = 1e-3
  stepsize = 1e-2
  iter = 0
  maxIter = 30
  while(iter < maxIter) {
    # browser()

    # fit the glm with the dispersions
    dispersionVector = numeric(length(groups))
    dispersionVector[groups == 1] = dispersions[1]
    dispersionVector[groups == 2] = dispersions[2]
    modelNegBin = suppressWarnings(glm.fit(modelMatrix, count,
                     #family = NBGeneral(dispersionVector),
                     family = MASS::negative.binomial( 1 / dispersionVector),
                     offset = log(countPerCell))
    )
    mu = modelNegBin$fitted.values
    logLik = -0.5 * (modelNegBin$aic -
                       2 * (length(count) - modelNegBin$df.residual))

    if (!is.null(logLikPrevious)) {
      if (logLik - logLikPrevious < eps ||
          max(abs(dispersionsPrevious - dispersions)) < epsDispersion) {
        break
      }
    }

    logLikPrevious = logLik
    dispersionsPrevious = dispersions

    # optimize the dispersions
    result = optim(dispersionsPrevious, likelihoodWithKnownMean, gr = NULL,
                   mu = mu, count = count, groups = groups,
                   modelMatrix = modelMatrix, type = "ML",
                   #method = "BFGS",
                   #lower = dispersionLowerBound, upper = 1e5,
                   control = list(fnscale = -1))
    dispersions = dispersionsPrevious +
                (result$par - dispersionsPrevious) * stepsize
    logLikOptim = result$value

    iter = iter + 1

    # <debug>
    # model = glm.fit(modelMatrix, count,
    #                       family = NBGeneral(dispersionVector),
    #                       #family = MASS::negative.binomial(dispersionVector),
    #                       offset = log(countPerCell))
    # model = glm(count ~ modelMatrix,
    #                 family = negative.binomial(1 / dispersionVector),
    #                 #family = MASS::negative.binomial(dispersionVector),
    #                 offset = log(countPerCell))
    # print(logLik(model))
    # # check with the direct calculation
    # logLikDirect = profileLogLikelihood2(dispersionVector, modelMatrix,
    #                                count, model$fitted.values, type = "ML")
    # print(logLikDirect)
    # </debug>

  }

  return(list(dispersions = dispersions, logLik = logLik, iter = iter))
}

MLETwoDispersions = function() {

}

likelihoodWithKnownDispersion = function() {
  # TODO
}

testCommonDispersion = function(count, groups,
                                              countPerCell,
                                              covariates,
                                              useGradient = F) {

  groups = factor(groups)

  # common mean structure
  if (!is.null(covariates)) {
    modelMatrix = model.matrix(~ groups + covariates)
  } else {
    modelMatrix = model.matrix(~ groups)
  }

  modelInit = suppressWarnings(glm.fit(modelMatrix, count,
                         family = MASS::negative.binomial( 1 / 0.1),
                         offset = log(countPerCell)))

  # one dispersion
  dispersionsGroupName1 = "oneGroup"
  commonGroups = rep(dispersionsGroupName1, length(count))
  otherData0 = list(count = count,
                   groups = commonGroups,
                   modelMatrix = modelMatrix,
                   countPerCell = countPerCell,
                   dispersionsGroupName = dispersionsGroupName1)
  parameterInitial0 = c(log(0.1), coefficients(modelInit))
  # get ML
  if (!useGradient) {
    gradientFun = NULL
    algorithm = "NLOPT_LN_NEWUOA"
  } else {
    gradientFun = gradientDifferentLogDispersionNegative
    algorithm = "NLOPT_LD_LBFGS"
  }
  opts = list(
    "algorithm" = algorithm,
    #"algorithm"="NLOPT_LN_NELDERMEAD",
    #"algorithm"="NLOPT_LN_NEWUOA",
    #"algorithm"="NLOPT_LN_BOBYQA",
    #"algorithm"="NLOPT_LN_COBYLA",
    #"algorithm"="NLOPT_LN_SBPLX",
    #"algorithm"="NLOPT_LD_MMA",
    #"algorithm"="NLOPT_LD_SLSQP",
    #"algorithm"="NLOPT_LD_LBFGS",
    #"algorithm"="NLOPT_LD_LBFGS_NOCEDAL",
    "xtol_rel"=1.0e-8,
    #"check_derivatives" = T,
    "maxeval" = 1e4)
   resultNLOPT0 = nloptr::nloptr(parameterInitial0,
                       eval_f = likelihoodDifferentLogDispersionNegative,
                       eval_grad_f = gradientFun,
                       opts = opts,
                       #lb = rep(-1e5, sum(!parameterPartition)),
                       #ub = rep(1e5, sum(!parameterPartition)),
                       otherData = otherData0)
   L0 = NA
   if (resultNLOPT0$status >= 0 && resultNLOPT0$status <= 4) {
     # MLEH0 = resultNLOPT0$solution
     L0 = -resultNLOPT0$objective
   } else {
     print("nloptr failed under one dispersion")
   }


  # # <debug>
  # # check gradiant function
  # browser()
  # gradNumeric = grad(likelihoodDifferentLogDispersionNegative,
  #                    MLEH0, otherData = otherData0)
  # gradTheoretic = gradientDifferentLogDispersionNegative(MLEH0,
  #                                             otherData = otherData0)
  #
  # # run with the gradient function
  # opts = list(
  #   #"algorithm" = algorithm,
  #   #"algorithm"="NLOPT_LN_NELDERMEAD",
  #   #"algorithm"="NLOPT_LN_NEWUOA",
  #   #"algorithm"="NLOPT_LN_BOBYQA",
  #   #"algorithm"="NLOPT_LN_COBYLA",
  #   #"algorithm"="NLOPT_LN_SBPLX",
  #   #"algorithm"="NLOPT_LD_MMA",
  #   #"algorithm"="NLOPT_LD_SLSQP",
  #   "algorithm"="NLOPT_LD_LBFGS",
  #   #"algorithm"="NLOPT_LD_LBFGS_NOCEDAL",
  #   "xtol_rel"=1.0e-8,
  #   #"check_derivatives" = T,
  #   "maxeval" = 1e4)
  # resultNLOPT0 = nloptr(parameterInitial0,
  #                       eval_f = likelihoodDifferentLogDispersionNegative,
  #                       eval_grad_f = gradientDifferentLogDispersionNegative,
  #                       opts = opts,
  #                       #lb = rep(-1e5, sum(!parameterPartition)),
  #                       #ub = rep(1e5, sum(!parameterPartition)),
  #                       otherData = otherData0)
  # # </debug>

  # two dispersion
  dispersionsGroupName2 = sort(unique(groups))
  nGroups = length(unique(groups))
  otherData1 = list(count = count,
                   groups = groups,
                   modelMatrix = modelMatrix,
                   countPerCell = countPerCell,
                   dispersionsGroupName = dispersionsGroupName2)
  parameterInitial1 = c(log(rep(0.1, nGroups)), coefficients(modelInit))
  # get ML
  resultNLOPT1 = nloptr:nloptr(parameterInitial1,
                        eval_f = likelihoodDifferentLogDispersionNegative,
                        eval_grad_f = gradientFun,
                        opts = opts,
                        #lb = rep(-1e5, sum(!parameterPartition)),
                        #ub = rep(1e5, sum(!parameterPartition)),
                        otherData = otherData1)
  MLEH1 = resultNLOPT1$solution
  L1 = NA
  if (resultNLOPT1$status >= 0 && resultNLOPT1$status <= 4) {
    L1 = -resultNLOPT1$objective
  } else {
    print("nloptr failed under two dispersions")
  }

  LRT = 2 * (L1 - L0)
  pvalueLog = NA
  if (!is.na(LRT)) {
    pvalueLog = pchisq(LRT, df = 1, lower.tail = F, log.p = T)
  }

  return(c(LRT, pvalueLog))
}

#' negative binomial models allowing independent dispersions
#'
#' This function tests the following hypothesis
#' H0: the same mu for different groups
#' H1: different mu for different groups
#' In both cases, the dispersion is group specific. The dispersion is estimated
#' only for the full model and then used for the reduced model
#'
#' @param count the vector of gene counts
#' @param groups the vector of group information
#' @param countPercell the total UMI of each cell
#' @param covariates the covariaites
#' @param sizeFactor the normalization factor, the effective count used will be
#'                   the countPerCell * sizeFactor. The size factor from
#'                   package scran needs to be divided by countPerCell first,
#'                   i.e., the size factor from scran is used as the effective
#'                   count directly.
#' @param singleDispersion whether to use one single dispersion for all cells.
#'        All groups will share this single dispersion parameter, i.e., no
#'        independent dispersions are used
#' @param dispMethod the method to estimate the dispersion

#' @return a vector of the following:
#'
#'          LR: likelihood rato
#'
#'          beta: the effect parameter
#'
#'          dispersionOutput: dispersions of different groups

#' @export
NBID = function(count, groups,
                countPerCell,
                covariates = NULL,
                sizeFactor = NULL,
                singleDispersion = F,
                dispMethod = "poisson-ML") {

  #browser()

  type = "ML"
  if (regexpr("CR", dispMethod) != -1) {
    type = "CR"
  }

  # H1:
  ## full model: group specific dispersions
  # estimate dispersion and calculate likelihood for each group

  # groups as a factor with sorting levels
  uniqueGroups = NULL
  if (!is.factor(groups)) {
    uniqueGroups = as.character(unique(sort(groups)))
    groups = factor(groups, levels = uniqueGroups)
  } else {
    uniqueGroups = levels(groups)
  }
  nGroup = length(uniqueGroups)
  dispersions = NULL
  logLikEachGroup = matrix(NA, nrow = nGroup)
  expectedMean = matrix(NA, nrow = nGroup)
  modelMatrix = NA

  LR = NA
  beta = NA
  dispersionOutput = rep(NA, nGroup)
  logLik0 = NA
  logLik1 = NA

  if (!is.null(covariates)) {
    covariatesAll = cbind(model.matrix(~groups), covariates)
  } else {
    covariatesAll = model.matrix(~groups)
  }
  covariatesAll = covariatesAll[, -1, drop = F]

  if (is.null(countPerCell)) {
    stop("countPerCell is missing")
  }
  if (!is.null(sizeFactor) && all(sizeFactor > 0)) {
    countPerCell = countPerCell * sizeFactor
  }
  #browser()
  result = NULL
  if (singleDispersion) {
    uniqueGroups = "single"
    result = MLENBDirect(count, rep(uniqueGroups, length(groups)),
                         covariates = covariatesAll,
                         countPerCell,
                         0.1, uniqueGroups, type)
  } else {
    result = MLENBDirect(count, groups,
                       covariates = covariatesAll,
                       countPerCell,
                       rep(0.1, length(uniqueGroups)), uniqueGroups, type)
  }

  beta = result$coef[2:nGroup]

  # # <debug>
  # print(result)
  # # </debug>
  #browser()

  dispersions = result$dispersions
  logLik1 = result$logLik

  if (any(is.na(dispersions))) {
    return(c(LR, beta, dispersionOutput))
  }

  # H0
  logLik0 = NA
  dispersions0 = c(NA, NA)
  # browser()
  test = "twodisp_fixed"
  if (test == "twodisp_fixed") {
    # use full model's estimated dispersion
    if (!is.null(covariates)) {
      modelMatrix = model.matrix(~covariates)
    } else {
      modelMatrix = matrix(1, length(groups), 1)
    }
    dispersionVector = numeric(length(groups))
    dispersionVector[ ] = dispersionLowerBound
    # browser()
    if (length(uniqueGroups) > 0) {
      if (length(uniqueGroups) == 1) {
        dispersionVector[ ] = dispersions
      } else {
        for (i in 1:length(uniqueGroups)) {
          dispersionVector[groups == uniqueGroups[i]] = dispersions[i]
        }
      }
    }
    try({
      modelH0 = suppressWarnings(glm.fit(modelMatrix, count,
                      family = MASS::negative.binomial(1 / dispersionVector),
                      offset = log(countPerCell)))
      #browser()
      logLik0 = -0.5 * (modelH0$aic - 2 * dim(modelMatrix)[2])
    })
  } else if (test == "twodisp_fullnull") {
    # This part might be broken

    groupCurrent = rep(1, length(count))
    modelFormula = count ~ groupCurrent
    modelFrame <- data.frame(groupCurrent)
    result0 = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                        countPerCell, "pooled-ML")
    dispersion0 = result0$dispersion
    logLik0 = result0$logLik
  } else {
    # resultH0 = MLEUnderH0EqualMean(count, groups, covariates, countPerCell,
    #                                dispersions)
    # browser()
    resultH0 = MLENBDirect(count, groups, covariates, countPerCell,
                           dispersions, uniqueGroups, type)
    logLik0 = resultH0$logLik
    dispersions0 = resultH0$dispersions
  }

  # likelihood ratio
  LR = 2 * (logLik1 - logLik0)

  #browser()
  dispersionOutput = numeric(length(uniqueGroups))
  dispersionOutput[] = dispersionLowerBound
  if (length(uniqueGroups) > 0) {
    dispersionOutput = dispersions
  }
  # <debug>
  #print("dispersions0")
  #print(dispersions0)
  # </debug>
  return(c(LR, beta, dispersionOutput))
}

NBIDNoCov = function(count, groups,
                countPerCell,
                covariates, dispMethod) {

  #browser()

  type = "ML"
  if (regexpr("CR", dispMethod) != -1) {
    type = "CR"
  }

  # H1:
  ## full model: group specific dispersions
  # estimate dispersion and calculate likelihood for each group
  uniqueGroups = as.character(unique(sort(groups)))
  nGroup = length(uniqueGroups)
  dispersions = NULL
  dispersionsGroupName = character(0)
  logLikEachGroup = matrix(NA, nrow = nGroup)
  expectedMean = matrix(NA, nrow = nGroup)
  modelMatrix = NA

  LR = NA
  beta = NA
  dispersionOutput = rep(NA, nGroup)
  logLik0 = NA
  logLik1 = NA

  for (i in 1:nGroup) {
    indexGroup = which(groups == uniqueGroups[i])
    groupCurrent = factor(groups[indexGroup])
    countCurrent = count[indexGroup]
    if (!is.null(covariates)) {
      covariatesCurrent = covariates[indexGroup, , drop = F]
      modelFormula = formula("count ~ covariatesCurrent")
      modelFrame <- list(groupCurrent = groupCurrent,
                         covariatesCurrent = covariatesCurrent)
      modelMatrix = model.matrix(~ covariatesCurrent)
    } else {
      modelFormula = formula("count ~ groupCurrent")
      modelFrame <- data.frame(groupCurrent)
      modelMatrix = matrix(1, length(countCurrent), 1)
    }

    result = MLENBDirect(countCurrent, groupCurrent,
                         covariates = covariates[indexGroup],
                         countPerCell[indexGroup],
                         c(0.1), uniqueGroups[i], type)
    expectedMean[i] = result$coef[1]

    # # <debug>
    # print(result)
    # # </debug>
    #browser()

    dispersions = c(dispersions, result$dispersions)
    dispersionsGroupName = c(dispersionsGroupName, uniqueGroups[i])
    logLikEachGroup[i] = result$logLik
  }
  logLik1 = sum(logLikEachGroup)

  if (any(is.na(dispersions))) {
    return(c(LR, beta, dispersionOutput, logLik0, logLik1))
  }

  # H0
  logLik0 = NA
  dispersions0 = c(NA, NA)
  # browser()
  test = "twodisp_fixed"
  if (test == "twodisp_fixed") {
    # use full model's estimated dispersion
    if (!is.null(covariates)) {
      modelMatrix = model.matrix(~covariates)
    } else {
      modelMatrix = matrix(1, length(groups), 1)
    }
    dispersionVector = numeric(length(groups))
    dispersionVector[ ] = dispersionLowerBound
    # browser()
    if (length(dispersionsGroupName) > 0) {
      for (i in 1:length(dispersionsGroupName)) {
        dispersionVector[groups == uniqueGroups[i]] = dispersions[i]
      }
    }
    try({
      modelH0 = suppressWarnings(glm.fit(modelMatrix, count,
                                         family = MASS::negative.binomial(1 / dispersionVector),
                                         offset = log(countPerCell)))
      logLik0 = -0.5 * (modelH0$aic - 2 * dim(modelMatrix)[2])
    })
  } else if (test == "twodisp_fullnull") {
    groupCurrent = rep(1, length(count))
    modelFormula = count ~ groupCurrent
    modelFrame <- data.frame(groupCurrent)
    result0 = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                        countPerCell, "pooled-ML")
    dispersion0 = result0$dispersion
    logLik0 = result0$logLik
  } else {
    # resultH0 = MLEUnderH0EqualMean(count, groups, covariates, countPerCell,
    #                                dispersions)
    # browser()
    resultH0 = MLENBDirect(count, groups, covariates, countPerCell,
                           dispersions, dispersionsGroupName, type)
    logLik0 = resultH0$logLik
    dispersions0 = resultH0$dispersions

    # # <debug>
    # browser()
    # beta = resultH0$coef
    # modelMatrix = matrix(1, length(count), 1)
    # mu = countPerCell * exp(modelMatrix %*% beta)
    # print("log likelihood for count value 0 and 1")
    # print(dnbinom(c(0, 1), mu = mu[1],
    #               size = 1 / dispersions0[1], log = T))
    # print("log likelihood for group 1")
    # logLikG1 = dnbinom(count[groups == 1], mu = mu[1],
    #                    size = 1 / dispersions0[1], log = T)
    # print(table(logLikG1))
    # print("sum of group1")
    # print(sum(logLikG1))
    # print("log likelihood for group 2")
    # logLikG2 = dnbinom(count[groups == 2], mu = mu[1],
    #                    size = 1 / dispersions0[2], log = T)
    # print(table(logLikG2))
    # print("sum of group2")
    # print(sum(logLikG2))
    # # </debug>
  }

  # likelihood ratio
  LR = 2 * (logLik1 - logLik0)

  #browser()
  beta = (expectedMean[2] - expectedMean[1])

  #browser()
  dispersionOutput = numeric(length(uniqueGroups))
  dispersionOutput[] = dispersionLowerBound
  if (length(dispersionsGroupName) > 0) {
    dispersionOutput[dispersionsGroupName %in% uniqueGroups] = dispersions
  }
  # <debug>
  #print("dispersions0")
  #print(dispersions0)
  # </debug>
  return(c(LR, beta, dispersionOutput))
}

#' Differential expression analysis using NBID on a data matrix
#'
#' @param data the gene-cell matrix
#' @param groups the group vector
#' @param covariates the covariates
#' @param countPerCell the total UMI per cell or the directly used reference
#'         count size. For example, this can be the size factor calculated
#'         from package scran
#' @param sizeFactor the normalization factor, the effective count used will be
#'                   the countPerCell * sizeFactor when both the countPerCell
#'                   and sizeFactor are set. In this case, the size factor from
#'                   package scran needs to be divided by countPerCell first
#'                   before assigned to this parameter to achieve the result
#'                   that the size factor from scran is used directly as the
#'                   count size.
#' @param singleDispersion whether to use one single dispersion for all cells.
#'        All groups will share this single dispersion parameter, i.e., no
#'        independent dispersions are used
#' @param dispMethod the method to estimate dispersions
#' @param ncore the number of cores to use for parallel running
#'
#' This function invokes NBID for each row of the data matrix.
#'
#'
#' @return a matrix of columns as follows:
#' pvalue: p value for each gene
#' LR: likelihood ratio test statistic
#' beta: the coefficient
#' dispersionGroup1, dispersionGroup2: estimated dispersions
#' log2FC: log2 fold change
#'
#' @examples
#' \dontrun{
#' data(smallData)
#' result = DEUsingNBID(smallData$count, smallData$groupLabel)
#' head(result)
#' }
#'
#' @export

DEUsingNBID = function(data, groups, covariates = NULL,
                          countPerCell = NULL,
                          sizeFactor = NULL,
                          singleDispersion = F,
                          dispMethod = "poisson-ML",
                          ncore = 1) {
  #browser()
  uniqueGroups = NULL
  if (!is.factor(groups)) {
    uniqueGroups = as.character(unique(sort(groups)))
    groups = factor(groups, levels = uniqueGroups)
  } else {
    uniqueGroups = levels(groups)
  }

  if (is.null(countPerCell)) {
    countPerCell = colSums(data)
  }
  if (is.null(covariates)) {
    design = model.matrix(~groups)
    design0 = model.matrix(~rep(1, dim(data)[2]) - 1)
  } else {
    design = model.matrix(~groups + covariates)
    design0 = model.matrix(~covariates)
  }

  stopifnot(length(unique(groups)) > 1)

  if (is.null(covariates)) {
    covariateString = "groups"
    covariateString0 = "1"
  } else {
    covariateString = "groups + covariates"
    covariateString0 = "covariates"
  }

  # start parallel running
  library(parallel)
  maxCores = detectCores()
  if (ncore > maxCores) {
    ncore = maxCores
    print(paste("maximum available cores is", maxCores))
  }
  cluster = makeCluster(ncore)
  print(paste("number of cores used for parallel computing", ncore))

  resultTest = parApply(cl = cluster, X = data, MARGIN = 1,
                     FUN = NBID, groups,
                     countPerCell, covariates, sizeFactor,
                     singleDispersion, dispMethod)


  stopCluster(cluster)
  # end parallel running

  nGroups = length(uniqueGroups)
  resultTest = t(resultTest)
  pvalue = pchisq(resultTest[, 1], df = nGroups - 1, lower.tail = F)
  log2FC = resultTest[, 2:nGroups] / log(2)
  if (!singleDispersion) {
    result = cbind(pvalue,
                 resultTest[, 1:(2 * nGroups)],
                 log2FC)
  } else {
    dispersions = matrix(NA, dim(resultTest)[1], nGroups)
    for (i in 1:nGroups) {
      dispersions[, i] = resultTest[, 3]
    }
    result = cbind(pvalue,
                   resultTest[, c(1, 2:nGroups)], dispersions,
                   log2FC)
  }
  colnames(result) = c("pvalue",
                       "LR",
                       paste0("beta", uniqueGroups[2:nGroups]),
                       paste0("dispersionGroup", uniqueGroups),
                       paste0("log2FC", uniqueGroups[2:nGroups]))
  return(result)
}
