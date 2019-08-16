#' @import pscl
#' @import MASS
#' @import nloptr
#' @import AER


strictLess = function(aic, delta) {
  # assume aic is ordered with model binomial, poisson, nb and zinb
  # choose the nb is only it is the minmum and is less than binomialby delta
  # choose the zinb model if it is the minimum and less than nb by delta

  minIndex = NA
  tryCatch({
    if (!is.na(aic[4]) && (min(aic[1:3], na.rm = T) - aic[4] > delta)) {
      minIndex = 4
    } else if (!is.na(aic[3]) && (min(aic[1:2]) - aic[3] > delta)) {
      minIndex = 3
    } else {
      minIndex = which.min(aic[1:2])
    }
  },
  error = function(e) {
    print(aic)
    print(e)
  }
  )


  return(minIndex)
}

checkCountContinuity = function(count) {
  # a outlier test to see whether there is an outlier
  # The idea is to table all the counts and see whether there is a
  # discontinuity in the table. The high count not connected to the lower count
  # is considered as outliers

  countTable = table(count)
  countValues = as.numeric(names(countTable))
  difference = diff(countValues)
  index = which(difference > 1)
  if (length(index) > 0 && index == length(difference) &&
      countTable[index + 1] == 1 &&
      countValues[index + 1] / countValues[index] > 2
      ) {
    print("raw")
    print(countTable)

    # set this count vaue to the closest smaller count
    indexCount = which(count == countValues[index + 1])
    count[indexCount] = countValues[index]

    print("modified")
    print(table(count))
    return(count)
  }
  return(count)
}

NBVsZINB = function(count, groupID, countPerCell,
                    covariates,
                    covariatesString = "1",
                    largestKRemoved = 0) {
    # fit different distribution to the count data
    # this function might not fit zinb very well,
    # so check the zinb logLik vs NB logLik after fitting to make sure
    # zinb is fitted well

    nCovariates = 0
    if (!is.null(covariates)) {
      covariates = as.matrix(covariates)
      p = dim(covariates)[2]
      if (p == 1 && length(unique(covariates)) == 1) {
        # only the intercept
        nCovariates = 1
      } else {
        nCovariates = p + 1
      }
    }

    if (largestKRemoved > 0) {
      indexRemove = order(count, decreasing = T)[1:largestKRemoved]
      count = count[-indexRemove]
      groupID = groupID[-indexRemove]
      countPerCell = countPerCell[-indexRemove]
      if (!is.null(covariates)) {
        covariates = covariates[-indexRemove, , drop = F]
      }
    }

    #browser()
    proportion = sum(count) / sum(countPerCell)
    maxCount = max(count)
    minCount = min(count)
    aveCount = mean(count)

    modelNB = list()
    modelZINB = list()

    modelNB$converged = NA
    modelZINB$converged = NA

    # theta estimated by the CR model
    aicCoxReid = NA
    dispersionCR = NA
    beta0CR = NA

    aicNB = NA
    thetaNB = NA
    thetaSENB = NA
    alphaNB = NA

    aicZINB = NA
    thetaZINB = NA
    logthetaSE = NA
    zeroProb = NA
    zeroPValue = NA

    diffAICNBZINB = NA
    diffLogLik = NA

    if (minCount > 0) {
      diffAICNBZINB = -Inf
      diffLogLik = 0
    } else {
      # negative binomial
      tic = proc.time()
      logLikNB = NA
      try({
        modelFormula = formula(paste("count ~", covariatesString))
        modelFrame <- data.frame(covariates)
        result0 = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                            countPerCell, "pooled-ML")
        dispersionCR = result0$dispersion
        logLikNB = result0$logLik
        # browser()
        beta0CR = coefficients(result0$model)[1]
        aicCoxReid = -2 * logLikNB + 2 * (nCovariates + 1)
        modelNB = result0$model
      })
      thetaNB = 1 / dispersionCR
      toc = proc.time()
      print("time of NBDispersionAndLikelihood")
      print(toc - tic)


      # #browser()
      # tic = proc.time()
      # initTheta = 1 #c(10^seq(from = -8, by = 1, to = 4)) # use 1 to save time
      # for (i in 1:length(initTheta)) {
      #   try({
      #     modelNBTry = glm.nb(formula(paste("count ~ ", covariatesString,
      #                                       "+ offset(log(countPerCell))")),
      #                         init.theta = initTheta[i],
      #                         control = glm.control(maxit = 30, epsilon = 1e-8))
      #     if (is.na(modelNB$converged) || logLik(modelNBTry) > logLik(modelNB)) {
      #       modelNB = modelNBTry
      #     }
      #   })
      # }
      #
      # # try with dispersionCR
      # if (!is.na(dispersionCR)) {
      #   try({
      #     modelNBTry = glm.nb(formula(paste("count ~ ", covariatesString,
      #                                       "+ offset(log(countPerCell))")),
      #                         init.theta = 1 / dispersionCR,
      #                         control = glm.control(maxit = 30, epsilon = 1e-8))
      #     if (is.na(modelNB$converged) || logLik(modelNBTry) > logLik(modelNB)) {
      #       modelNB = modelNBTry
      #     }
      #   })
      # }
      # if (!is.na(modelNB$converged)) {
      #   aicNB = modelNB$aic
      #   thetaNB = modelNB$theta
      #   thetaSENB = modelNB$SE.theta
      #   alphaNB = 1 / modelNB$theta
      # }
      #
      # toc = proc.time()
      # print("time of glm.nb")
      # print(toc - tic)

      # zinb
      tic = proc.time()
      try({
        # modelZINB = zeroinfl(formula(paste("count ~ 1 | 1")),
        #                         dist = "negbin", offset = log(countPerCell),
        #                   EM = T,  maxit = 500)
        #browser()
        if (!is.na(dispersionCR)) {
          EM = F
          start = list(count = coefficients(modelNB), zero = -3, theta = thetaNB)
        } else {
          EM = F
          start = NULL
        }
        modelZINB = pscl::zeroinfl(
                       formula(paste("count ~", covariatesString, "| 1")),
                       dist = "negbin", offset = log(countPerCell),
                       EM = EM,
                       start = start)
      })
      if (is.na(modelZINB$converged)) {
        try({
          modelZINB = pscl::zeroinfl(
                        formula(paste("count ~", covariatesString, "| 1")),
                        dist = "negbin", offset = log(countPerCell),
                        EM = T)
        })
      }
      toc = proc.time()
      print("time of zinb")
      print(toc - tic)

      logLikZINB = NA
      if (!is.na(modelZINB$converged)) {
        logLikZINB = logLik(modelZINB)
        aicZINB = -2 * logLikZINB + 2 * (nCovariates + 2)
        thetaZINB = modelZINB$theta
        logthetaSE = modelZINB$SE.logtheta
        zeroProb = 1 / (1 + exp(-modelZINB$coefficients$zero))
        zeroPValue = summary(modelZINB)$coefficients$zero[4]
      }

      # if (aicNB < aicCoxReid - 0.1) {
      #   print(c(aicNB, aicCoxReid))
      #   browser()
      # }
      # diffAICNBZINB = min(aicNB, aicCoxReid) - aicZINB
      diffAICNBZINB = aicCoxReid - aicZINB
      diffLogLik = logLikZINB - logLikNB
    }
    #output
    # browser()
    result = data.frame(proportion, minCount, maxCount, aveCount,
                        aicNB,
                        aicCoxReid,
                        diffAICNBZINB,
                        diffLogLik,
                        dispersionCR, alphaNB, thetaNB, thetaSENB,
                        thetaZINB, logthetaSE, zeroProb, zeroPValue,
                        modelNB$converged,
                        modelZINB$converged)

    return(result)
}

#' fit different count models
#'
#' This function fits the count with different count models including Poisson,
#' negative binomial and zero inflated negative binomial models
#'
#' @param count the vector of counts
#' @param groupID the group ID character
#' @param countPerCell the total UMI per cell
#' @param covariates the covariates to account for
#' @param covariatesString the covariates string used in the formula
#' @param largestKRemoved remove the numbers correspond to the largest K
#'
#' @return a data frame with fitted information
#'
#' @examples
#' \dontrun{
#' count = rnbinom(100, size = 2, mu = 0.1)
#' groupID = 1
#' countPerCell = floor(rnorm(100, 1e4, 3000))
#' fitResult = checkCountDistributionPerGroup(count, groupID, countPerCell)
#' print(t(fitResult))
#' }
#'
#' @export
checkCountDistributionPerGroup = function(count, groupID, countPerCell,
                                          covariates = NULL,
                                          covariatesString = "1",
                                          largestKRemoved = 0) {
  # fit different distribution to the count data

  if (is.null(covariates) || unique(covariates) == 1) {
    covariatesString = "1"
    covariates = rep(1, length(count))
  } else {
    covariatesString = "covariates"
  }

  nCovariates = 0
  if (!is.null(covariates)) {
    covariates = as.matrix(covariates)
    p = dim(covariates)[2]
    if (p == 1 && length(unique(covariates)) == 1) {
      # only the intercept
      nCovariates = 1
    } else {
      nCovariates = p + 1
    }
  }

  if (largestKRemoved > 0) {
    indexRemove = order(count, decreasing = T)[1:largestKRemoved]
    count = count[-indexRemove]
    groupID = groupID[-indexRemove]
    countPerCell = countPerCell[-indexRemove]
    if (!is.null(covariates)) {
      covariates = covariates[-indexRemove, , drop = F]
    }
  }

  #browser()
  proportion = sum(count) / sum(countPerCell)
  maxCount = max(count)
  minCount = min(count)
  aveCount = mean(count)

  modelBinomial = list()
  modelPoisson = list()
  modelQB = list()
  modelQP = list()
  modelNB = list()
  modelZINB = list()

  modelBinomial$converged = NA
  modelPoisson$converged = NA
  modelQB$converged = NA
  modelQP$converged = NA
  modelNB$converged = NA
  modelZINB$converged = NA

  pvalueDispLRT = NA

  # binomial
  aicBinomial = NA
  try({
  modelBinomial = glm(formula(paste("count / countPerCell ~ ",
                                    covariatesString)),
      weights = countPerCell,
      family = "binomial", control = glm.control(maxit = 30, epsilon = 1e-8))
  # modelBinomial = glm(formula(paste("cbind(count, countPerCell - count) ~
  #                          covariatesString")),
  #    family = "binomial", control = glm.control(maxit = 30, epsilon = 1e-8))
  aicBinomial = modelBinomial$aic
  })

  # poisson
  aicPoisson = NA
  try({
  modelPoisson = glm(formula(paste("count ~ ", covariatesString)),
                     offset = log(countPerCell),
      family = "poisson",
      control = glm.control(maxit = 30, epsilon = 1e-8)
      )
  aicPoisson = modelPoisson$aic
  })

  # quasi binomial
  aicQB = NA
  try({
  modelQB = glm(formula(paste("count / countPerCell ~ ", covariatesString)),
                      weights = countPerCell,
  family = "quasibinomial", control = glm.control(maxit = 30, epsilon = 1e-8))
  dispQB = summary(modelQB)$dispersion
  aicQB = modelQB$aic
  })

  # quasi poisson
  aicQP = NA
  try({
  modelQP = glm(formula(paste("count ~ ", covariatesString)),
                     offset = log(countPerCell),
    family = "quasipoisson", control = glm.control(maxit = 30, epsilon = 1e-8))
  dispQP = summary(modelQP)$dispersion
  aicQP = modelQP$aic
  })

  #### negative binomial
  # theta estimated by the CR model
  aicCoxReid = NA
  dispersionCR = NA
  beta0CR = NA
  modelNBFullCR = list()
  modelNBFullCR$converged = NA
  try({
    modelFormula = formula(paste("count ~", covariatesString))
    modelFrame <- data.frame(covariates)
    modelNBFullCR = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                        countPerCell, "full-CR")
    if (!modelNBFullCR$converged) {
      print("not converged")
    }
    dispersionCR = modelNBFullCR$dispersion
    logLikCR = modelNBFullCR$logLik
    # browser()
    beta0CR = coefficients(modelNBFullCR$model)[1]
    aicCoxReid = -2 * logLikCR + 2 * (nCovariates + 1)
  })

  #browser()
  aicPooledML = NA
  dispersionPooledML = NA
  beta0PooledML = NA
  modelNBPooledML = list()
  modelNBPooledML$converged = NA
  try({
    modelFormula = formula(paste("count ~", covariatesString))
    modelFrame <- data.frame(covariates)
    modelNBPooledML = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                        countPerCell, "pooled-ML")
    if (!modelNBPooledML$converged) {
      print("not converged")
    }
    dispersionPooledML = modelNBPooledML$dispersion
    logLikPooledML = modelNBPooledML$logLik
    # browser()
    beta0PooledML = coefficients(modelNBPooledML$model)[1]
    aicPooledML = -2 * logLikPooledML + 2 * (nCovariates + 1)
  })

  aicNB = NA
  thetaNB = NA
  thetaSENB = NA
  alphaNB = NA
  #browser()
  initTheta = c(10^seq(from = -8, by = 1, to = 4))
  for (i in 1:length(initTheta)) {
    try({
      modelNBTry = suppressWarnings(
            MASS::glm.nb(formula(paste("count ~ ", covariatesString,
                                    "+ offset(log(countPerCell))")),
                          init.theta = initTheta[i],
                     control = glm.control(maxit = 30, epsilon = 1e-8))
      )
      if (is.na(modelNB$converged) || logLik(modelNBTry) > logLik(modelNB)) {
        modelNB = modelNBTry
      }
    })
  }
  #browser()
  # try with dispersionCR
  if (!is.na(dispersionCR)) {
    try({
      modelNBTry = suppressWarnings(
                MASS::glm.nb(formula(paste("count ~ ", covariatesString,
                                        "+ offset(log(countPerCell))")),
                          init.theta = 1 / dispersionCR,
                          control = glm.control(maxit = 30, epsilon = 1e-8))
      )
      if (is.na(modelNB$converged) || logLik(modelNBTry) > logLik(modelNB)) {
        modelNB = modelNBTry
      }
    })
  }
  if (!is.na(modelNB$converged)) {
    aicNB = modelNB$aic
    thetaNB = modelNB$theta
    thetaSENB = modelNB$SE.theta
    alphaNB = 1 / modelNB$theta
  }

  # zinb
  # set EM to T has the risk of non-stop running in zeroinfl!
  aicZINB = NA
  thetaZINB = NA
  logthetaSE = NA
  zeroProb = NA
  zeroPValue = NA
  if (minCount > 0) {
    aicZINB = aicNB + 2
    thetaZINB = NA
    logthetaSE = NA
    zeroProb = 0
    zeroPValue = NA
    modelZINB$converged = T
  } else {
    try({
    # modelZINB = zeroinfl(formula(paste("count ~ 1 | 1")),
    #                         dist = "negbin", offset = log(countPerCell),
    #                   EM = T,  maxit = 500)
      #browser()
      if (!is.na(thetaNB)) {
        EM = F
        start = list(count = coefficients(modelNB), zero = -3, theta = thetaNB)
      } else {
        EM = F
        start = NULL
      }
      modelZINB = suppressWarnings(
                   pscl::zeroinfl(
                         formula(paste("count ~", covariatesString, "| 1")),
                         dist = "negbin", offset = log(countPerCell),
                         EM = EM,
                         start = start)
      )
    })
    if (is.na(modelZINB$converged)) {
      try({
        modelZINB = suppressWarnings(
            pscl::zeroinfl(
                             formula(paste("count ~", covariatesString, "| 1")),
                             dist = "negbin", offset = log(countPerCell),
                             EM = F)
        )
      })
    }
    if (!is.na(modelZINB$converged)) {
      aicZINB = -2 * logLik(modelZINB) + 2 * (nCovariates + 2)
      thetaZINB = modelZINB$theta
      logthetaSE = modelZINB$SE.logtheta
      zeroProb = 1 / (1 + exp(-modelZINB$coefficients$zero))
      zeroPValue = summary(modelZINB)$coefficients$zero[4]
    }
  }

  # # fit the Poisson tweedie model
  # browser()
  # countMatrix = matrix(count, nrow = 1)
  # rownames(countMatrix) = "gene"
  # aicPT = NA
  # modelPTConverged = F
  # muPT = NA
  # dPT = NA
  # aPT = NA
  # try({
  #   resultPT = PTmodel(countMatrix, 0, countPerCell)$StatMatrix
  #   aicPT = resultPT[, "AIC"]
  #   modelPTConverged = is.null(resultPT[, "Convergence"])
  #   muPT = resultPT[, "mu"]
  #   dPT = resultPT[, "D"]
  #   aPT = resultPT[, "a"]
  # })

  diffAICNBZINB = aicNB - aicZINB
  diffAICBinomialNB = aicBinomial - aicNB

  # best model based on aic
  aicNBBest = apply(cbind(aicPooledML, aicNB), 1, min, na.rm = T)
  aics = cbind(aicBinomial, aicPoisson, aicNBBest, aicZINB)
  modelNames = c("binomial", "poisson", "nb", "zinb")
  modelByAIC = modelNames[which.min(aics)]

  # poisson vs NB test based on dispersiontest from AER
  pvalueDispTest = AER::dispersiontest(modelPoisson)$p.value

  # poisson vs NB test overdispersion test based on LRT
  # assign("count", count, envir = .GlobalEnv)
  # assign("groups", groups, envir = .GlobalEnv)
  # assign("countPerCell", countPerCell, envir = .GlobalEnv)
  # assign("covariateString", covariateString, envir = .GlobalEnv)
  # pvalueodTest = odTest(modelNB)
  try({
    LRTScore = 2*(logLik(modelNB) - logLik(modelPoisson))
    if (LRTScore <= 0) {
      pvalueDispLRT = 1
    } else {
      pvalueDispLRT = 0.5 * pchisq(LRTScore, df = 1, lower.tail = F)
    }
  })

  #output
  # browser()
  result = data.frame(proportion, minCount, maxCount, aveCount,
                 aicCoxReid, aicPooledML,
                 aicBinomial, aicPoisson, aicNB, aicZINB,
                 #aicPT,
                 modelByAIC,
                 diffAICNBZINB, diffAICBinomialNB,
                 dispQB, dispQP,
                 dispersionCR, alphaNB, thetaNB, thetaSENB,
                 thetaZINB, logthetaSE, zeroProb, zeroPValue,
                 pvalueDispTest, pvalueDispLRT, beta0CR,
                 #muPT, dPT, aPT,
                 #modelPTConverged,
                 modelBinomial$converged, modelPoisson$converged,
                 modelQB$converged, modelQP$converged, modelNB$converged,
                 modelZINB$converged,
                 modelNBFullCR$converged,
                 modelNBPooledML$converged
                 )

  # add group ID
  colnames(result) = paste0(groupID, "_", colnames(result))
  rownames(result) = NULL

  # #<debug>
  # #browser()
  # print(t(result))
  # #</debug>

  return(result)
}

checkCountDistribution = function(count, groups, countPerCell,
                                  covariates = NULL, covariatesString = NULL,
                                  largestKRemoved = 0) {
  # check the distribution in each group

  result = NULL
  uniqueGroups = unique(groups)
  for (i in 1:length(uniqueGroups)) {
    groupID = uniqueGroups[i]
    indexPerGroup = which(groups == groupID)
    countPerGroup = count[indexPerGroup]
    countPerCellPerGroup = countPerCell[indexPerGroup]
    resultPerGroup = checkCountDistributionPerGroup(countPerGroup, groupID,
                                                    countPerCellPerGroup,
                                                    covariates,
                                                    covariatesString,
                                                    largestKRemoved)
    if (is.null(result)) {
      result = resultPerGroup
    } else {
      result = cbind(result, resultPerGroup)
    }
  }

  #browser()
  return(result)
}

#' fit different count models for a data matrix
#'
#' This function fits the count with different count models including Poisson,
#' negative binomial and zero inflated negative binomial models. It invokes
#' checkCountDistributionPerGroup for the real fitting
#'
#' @param data the gene-cell data matrix, it assumes the row names of the matrix
#'         are the gene names
#' @param groups the group ID vector. If there are multiple groups, the results
#'        will be put together with cbind
#' @param countPerCell the total UMI per cell. If it is NULL, it is calculated
#'         by summing over all the gene counts for each cell.
#' @param covariates the covariates to account for
#' @param largestKRemoved remove the numbers correspond to the largest K
#' @param ncore the number of cores to use for parallel running
#'
#' @return a data frame with fitted information
#'
#' @examples
#' \dontrun{
#' gene1 = rnbinom(100, size = 0.1, mu = 0.1)
#' gene2 = rpois(100, lambda = 0.1)
#' data = rbind(gene1, gene2)
#' groups = rep(1, 100)
#' countPerCell = rep(1e4, 100)
#' fitResult = checkDataDistribution(data, groups, countPerCell)
#' print(t(fitResult))
#' }
#'
#' @export
checkDataDistribution = function(data, groups, countPerCell = NULL,
                                 covariates = NULL,
                                 largestKRemoved = 0,
                                 ncore = 1) {
  #browser()

  if (is.null(covariates)) {
    covariatesString = "1"
    covariates = rep(1, dim(data)[2])
  } else {
    covariatesString = "covariates"
  }

  if (is.null(countPerCell)) {
    countPerCell = colSums(data)
  }

  # keep rows with at least one nonzero count
  minCount = 1
  indexKept = (rowSums(data) >= minCount)
  data = data[indexKept, , drop = F]

  # result = apply(X = data, MARGIN = 1, FUN = checkCountDistribution,
  #                groups,
  #                countPerCell)

  # start parallel running
  library(parallel)
  maxCores = detectCores()
  if (ncore > maxCores) {
    ncore = maxCores
    print(paste("maximum available cores is", maxCores))
  }
  cluster = makeCluster(ncore)
  print(paste("number of cores used for parallel computing", ncore))

  result = parLapply(cl = cluster,
                  X = data.frame(t(data)), fun = checkCountDistribution,
                 groups,
                 countPerCell, covariates, covariatesString, largestKRemoved)

  stopCluster(cluster)
  # end parallel running

  result = do.call(rbind.data.frame, result)

  geneNames = rownames(data)
  result = cbind(geneNames, result)

  return(result)
}

#' summarize the model comparison results
#'
#' This function summarizes the results and count the number of models that are
#' selected for each gene. The selection criterion can be hypothesis testing
#' based, or AIC based.
#'
#' @param result
#'
#' @return a list of two components:
#'      modelCount: summarized count of best fitted models
#'      modelStats: the FDR and p-values of hypothesis testing results
#'
#' @examples
#' \dontrun{
#' gene1 = rnbinom(100, size = 0.1, mu = 0.1)
#' gene2 = rpois(100, lambda = 0.1)
#' data = rbind(gene1, gene2)
#' groups = rep(1, 100)
#' countPerCell = rep(1e4, 100)
#' fitResult = checkDataDistribution(data, groups, countPerCell)
#' summaryResult = summaryComparison(fitResult)
#' print(summaryResult)
#' }
#'
#' @export
summaryComparison = function(result) {
  # change row names to remove the group label
  colnames(result) = sub("^.*_", "", colnames(result))

  # check the number of genes with converged poisson, NB and ZINB
  indexConverged = which(
    result[["modelPoisson.converged"]] &
      (result[["modelNB.converged"]] | result[["modelNBPooledML.converged"]]) &
      (result[["modelZINB.converged"]])
  )
  nConverged = length(indexConverged)
  nExamined = dim(result)[1]
  #print(paste(nConverged,
  #            "number of genes converged out of", nExamined))

  # analyze only the converged ones
  result = result[indexConverged, , drop = F]

  threshold = 0.5 # numerical precision threshold when comparing two logLik

  # make sure logLik of zinb is not smaller than nb
  # so no obvious fitting issues: check that Max > -1
  aicZINB = result[["aicZINB"]]
  aicPooledML = result[["aicPooledML"]]
  aicNB = result[["aicNB"]]
  aicNBUsed = apply(cbind(aicNB, aicPooledML), 1, min, na.rm = T)
  aicBinomial = result[["aicBinomial"]]
  aicPoisson = result[["aicPoisson"]]
  LRTDifference = -(aicZINB - (aicNBUsed + 2)) / 2
  #stopifnot(min(LRTDifference, na.rm = T) > -1)
  indexToCheck = which(LRTDifference <= -1 * threshold)

  # # debug information
  # if (length(indexToCheck) > 0) {
  #   print(paste(length(indexToCheck),
  #               "genes whose logLik of zinb is", threshold, "smaller than NB.",
  #               "This shows problems when fitting zinb"))
  #   if (length(indexToCheck) > 10) {
  #     print("at most 10 genes printed")
  #   }
  #   print(result[indexToCheck[1:min(10, length(indexToCheck))], ])
  # }
  # ####

  # check min(logLik(NB) - logLik(Binomial)), for most cases, NB should
  # have a better fit than binomial, but not gauranteed as compared to
  # poisson
  # print("summary(logLik(NBUsed) - logLik(Binomial))")
  # print(summary(-(aicNBUsed - (aicBinomial + 2)) / 2, na.rm = T))
  # if (min(-(aicNBUsed - (aicBinomial + 2)) / 2, na.rm = T) <= -1 * threshold) {
  #   index = which(-(aicNBUsed - (aicBinomial + 2)) / 2 <= -1 * threshold)
  #   print(paste("# of genes (binom > NB):", length(index)))
  #   #indexToCheck = union(indexToCheck, index)
  # }

  # check min(logLik(NB) - logLik(poisson)) > -0.5
  # print(paste("min(logLik(NBUsed) - logLik(Poisson))",
  #             min(-(aicNBUsed - (aicPoisson + 2)) / 2, na.rm = T)))
  # print(sum(-(aicNBUsed - (aicPoisson + 2)) / 2 < -1 * threshold))
  # #stopifnot(min(-(aicNBUsed - (aicPoisson + 2)) / 2, na.rm = T) >
  # #            -1 * threshold)

  LRTDifference = -(aicNBUsed - (aicPoisson + 2)) / 2
  indexToCheck2 = which(LRTDifference <= -1 * threshold)

  # if (length(indexToCheck2) > 0) {
  #   print(paste(length(indexToCheck2),
  #               "genes whose logLik of NB is", threshold,
  #               "smaller than poisson.",
  #               "This shows problems when fitting NB"))
  #   if (length(indexToCheck2) > 10) {
  #     print("at most 10 genes printed")
  #   }
  #   print(result[indexToCheck2[1:min(10, length(indexToCheck2))], ])
  # }

  indexToCheck = c(indexToCheck, indexToCheck2)
  if (length(indexToCheck) > 0) {
    result = result[-indexToCheck, , drop = F]
  }
  nUsed = dim(result)[1]
  print(paste("genes passing convergence check:", nUsed))

  minCount = result[["minCount"]]
  geneNames = result[[1]]
  bestModel = result[["modelByAIC"]]
  aicNB = result[["aicNB"]]
  aicBinomial = result[["aicBinomial"]]
  aicPoisson = result[["aicPoisson"]]
  aicZINB = result[["aicZINB"]]
  aicCoxReid = result[["aicCoxReid"]]
  aicPooledML = result[["aicPooledML"]]

  aicNBUsed = apply(cbind(aicNB, aicPooledML), 1, min, na.rm = T)

  stopifnot(sum(is.na(aicBinomial)) == 0)
  stopifnot(sum(is.na(aicPoisson)) == 0)
  #stopifnot(sum(is.na(aicNB)) == 0)
  #stopifnot(sum(is.na(aicZINB[minCount == 0])) == 0)
  #print("no missing aics in binomial, poisson, NB and ZINB")

  # check logLik(poisson) - logLik(binomial) < 1 when the best model is
  # poisson or binomial
  #stopifnot(abs(aicPoisson - aicBinomial)[
  #             bestModel == "poisson" | bestModel == "binomial"] / 2 < 1)

  # print("logLik poisson - binomial when the best model is poisson or binom")
  # print(summary((aicPoisson - aicBinomial)[
  #   bestModel == "poisson" | bestModel == "binomial"] / (-2)))
  #
  # print("Ideally aic(NBPooledML) - aic(NB) is not much larger than 0")
  # print(summary(aicPooledML - aicNB))
  #
  # print("Ideally aic(NB) - aic(NBCR) is not much larger than 0")
  # print(summary(aicNB - aicCoxReid))
  #
  # print("Ideally aic(NBPooledML) - aic(NBCR) is not much larger than 0")
  # print(summary(aicPooledML - aicCoxReid))
  #
  # print("make sure aic(NBUsed) - aic(NBCR) is not much larger than 0")
  # print(summary(aicNBUsed - aicCoxReid))

  aic = cbind(aicPoisson, aicNBUsed, aicZINB)
  delta = 0
  bestModelStrictIndex = apply(X = aic, MARGIN = 1, FUN = strictLess, delta)
  modelNames = c("poisson", "nb", "zinb")
  bestModelStrict = modelNames[bestModelStrictIndex]
  # print("best model table")
  # print(table(bestModel))
  # print("best model table with AIC")
  # print(table(bestModelStrict))

  # browser()

  diffNBZINB = aicNBUsed - aicZINB

  # do a likelihood ratio test between NB and ZINB
  LRTStat = (aicNBUsed + 2) - aicZINB
  pvalueNBVsZINB = 0.5 * pchisq(LRTStat, df = 1, lower.tail = F)
  fdrNBVsZINB = p.adjust(pvalueNBVsZINB, method = "BH")

  #browser()
  # do a likelihood ratio test between poisson and NB
  LRTStat = (aicPoisson + 2) - aicNBUsed
  pvaluePoissonVsNB = 0.5 * pchisq(LRTStat, df = 1, lower.tail = F)
  pvaluePoissonVsNBNoZINB = pvaluePoissonVsNB
  pvaluePoissonVsNBNoZINB[fdrNBVsZINB < 0.05] = NA
  fdrPoissonVsNBNoZINB = p.adjust(pvaluePoissonVsNBNoZINB, method = "BH")
  # assign models based on the stepwise test
  modelBasedOnTest = character(length(aicNBUsed))
  modelBasedOnTest[ ] = "poisson"
  modelBasedOnTest[fdrNBVsZINB < 0.05] = "zinb"
  modelBasedOnTest[!is.na(fdrPoissonVsNBNoZINB) &
                     fdrPoissonVsNBNoZINB < 0.05] = "nb"
  print(paste("model based on hypothesis test"))
  print(table(modelBasedOnTest))

  # print(paste("number of ZINB with FDR < 0.05:",
  #             sum(fdrNBVsZINB < 0.05, na.rm = T)))

  output = data.frame(fdrNBVsZINB, pvalueNBVsZINB, LRTStat,
                 diffNBZINB, bestModelStrict,
                 pvaluePoissonVsNB, fdrPoissonVsNBNoZINB,
                 modelBasedOnTest)
  modelStats = output[order(output[, 2]), , drop = F]

  # summary output
  modelCount = cbind(nExamined, nUsed,
                     sum(fdrNBVsZINB < 0.05, na.rm = T),
                     sum(modelBasedOnTest == "nb"),
                     sum(modelBasedOnTest == "poisson"),
          sum(modelBasedOnTest == "poisson") / length(modelBasedOnTest) * 100,
                     sum(bestModelStrict == "poisson"),
                     sum(bestModelStrict == "nb"),
                     sum(bestModelStrict == "zinb")
  )
  colnames(modelCount) = c("#(genes)", "converged", "test_zinb",
                           "test_NB", "test_Poisson", "percentage_possion",
                           "aic_Poisson", "aic_NB", "aic_ZINB"
  )

  return(list(modelCount = modelCount, modelStats = modelStats))
}
