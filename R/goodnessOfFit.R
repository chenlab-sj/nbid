#' @importFrom fitdistrplus plotdist
#' @importFrom osDesign rmvhyper

downSampleVector = function(vector, newTotalCounts) {
  # down sample a vector with a new total counts
  # if the new total counts is the same as sum(vector), then no down sampling
  # is performed

  count = numeric(length(vector))
  count[ ] = NA
  if (sum(!is.na(vector)) > 0) {
    count = rmvhyper(vector, newTotalCounts)
  }
  return(count)
}

#' downsample a data matrix
#'
#' @param data the gene*cell count matrix
#' @param groupLabel the label of each cell
#' @param newTotalCounts the UMI count to after downsampling
#' @param newQuantile the common UMI after downsampling is set by the quantile
#'             of the total UMI across cells before downsampling
#'             If newTotalCounts is set, this parameter will be ignored
#'
#' @return A list of two components
#'         data: this is the count matrix
#'         groupLable: this is the group label of cells
#'
#' @export
downsampleData = function(data, groupLabel,
                         newTotalCounts = NULL,
                         newQuantile = 0.25) {
  # for downsample, keep those cells with count >= newTotalCounts

  # input
  # data: a gene by cell matrix

    countPerCell = colSums(data)
    if (is.null(newTotalCounts) && !(is.null(newQuantile))) {
      newTotalCounts = quantile(countPerCell, probs = newQuantile)
    }
    cellRemove = which(countPerCell < newTotalCounts)
    if (length(cellRemove) > 0 ) {
      data = data[, -cellRemove]
      groupLabel = groupLabel[-cellRemove]
    }
    geneNames = rownames(data)
    data = apply(X = data, MARGIN = 2, FUN = downSampleVector, newTotalCounts)
    rownames(data) = geneNames


  return(list(data = data, groupLabel = groupLabel))
}


pzinb = function(x, mu, size, zeroProb) {
  # calculate the probability of zero inflated negative binomial distribution

  # input
  # x: a vector to evaluate the probability
  # mu: the mu of the NB component
  # size: the size of the NB component, it is 1 / dispersion
  # zeroProb: the probability of the zero component

  prob = NA
  prob = zeroProb +
    (1 - zeroProb) * pnbinom(x, size = size, mu = mu)
  prob[x < 0] = 0
  return(prob)
}

fitZINB = function(count) {
  # use NB fit as an initial fit and give it to zeroinfl for further fitting
  # assume no offset for all counts

  # input
  # count: a vector of counts

  # output
  # a list of three components:
  #   fitted: the fitted probablity
  #   diffLogLikZINBNB: the difference of loglik between zinb and nb
  #   converge: if succeed, it is "converged", else other messages
  #   modelZINB: the zero inflated model

  fitted = NA
  diffLogLikZINBNB = NA
  converge = "failed"
  modelZINB = NULL
  zeroProb = NA
  # maker sure there is at least one 0
  if (min(count) != 0) {
    converge = "no zero count"
    return(list(fitted = fitted, diffLogLikZINBNB = diffLogLikZINBNB,
                converge = converge, modelZINB = modelZINB,
                zeroProb = zeroProb))
  }

  modelNB = list()
  modelZINB = list()

  modelNB$converged = NA
  modelZINB$converged = NA

  # fit with NB
  logLikNB = NA
  countPerCell = numeric(length(count))
  countPerCell[ ] = 1
  try({
    modelFormula = formula(paste("count ~ 1"))
    modelFrame = data.frame(rep(1, length(count)))
    result0 = NBDispersionAndLikelihood(count, modelFormula, modelFrame,
                                        countPerCell, "pooled-ML")
    dispersionCR = result0$dispersion
    logLikNB = result0$logLik
    # browser()
    beta0CR = coefficients(result0$model)[1]
    modelNB = result0$model
  })
  thetaNB = 1 / dispersionCR

  try({
    if (!is.na(dispersionCR)) {
      EM = F
      start = list(count = coefficients(modelNB), zero = -3, theta = thetaNB)
    } else {
      EM = F
      start = NULL
    }
    modelZINB = suppressWarnings(
                   zeroinfl(formula(paste("count ~ 1 | 1")),
                         dist = "negbin",
                         EM = EM,
                         start = start)
    )
  })
  if (is.na(modelZINB$converged)) {
    try({
      modelZINB = suppressWarnings(
          zeroinfl(formula(paste("count ~1 | 1")),
                           dist = "negbin",
                           EM = T)
      )
    })
  }

  diffLogLikZINBNB = logLik(modelZINB) - logLikNB
  if (!is.na(modelZINB$converged) && !is.na(diffLogLikZINBNB) &&
      diffLogLikZINBNB > -0.5) {
    converge = "converged"
    thetaZINB = modelZINB$theta
    logthetaSE = modelZINB$SE.logtheta
    zeroProb = 1 / (1 + exp(-modelZINB$coefficients$zero))
    zeroPValue = summary(modelZINB)$coefficients$zero[4]
  }

  return(list(fitted = fitted, diffLogLikZINBNB = diffLogLikZINBNB,
              converge = converge, modelZINB = modelZINB,
              zeroProb = zeroProb))
}

#' goodness of fit for a count vector
#'
#' The distribution can be Poisson, negetive binomial or zero inflated
#' negative binomial (ZINB). For Poisson, glm is used for the fitting
#' For NB, glm.nb is used for the fitting. For ZINB,
#' a NB is first fitted and then zeroinfl is used for fitting
#'
#' @param count a count vector
#' @param distribution "poisson" for Poisson distribution,
#'        "nb" for the negative binomial distribution, "zinb" for the zero
#'        inflated negative binomial distribution
#' @param plot If TRUE, plot the empirical and theoretical density and
#'        cumulative distribution function using package fitdistrplus
#'
#' @return a data frame with the pvalue of the goodness of fit and whether the
#'         fit is converged
#'
#' @examples
#' \dontrun{
#' gene1 = rnbinom(500, size = 0.1, mu = 1)
#' print(goodnessOfFit(gene1, "poisson"))
#' print(goodnessOfFit(gene1, "nb"))
#' print(goodnessOfFit(gene1, "zinb"))
#'
#' gene2 = rpois(500, lambda = 1)
#' print(goodnessOfFit(gene2, "poisson"))
#' print(goodnessOfFit(gene2, "nb"))
#' print(goodnessOfFit(gene2, "zinb"))
#' }
#'
#' @export

goodnessOfFit = function(count, distribution, plot = F) {

  # fit the count data with the distribution
  fitResult = NULL
  data = data.frame(count)
  size = NULL
  mu = NULL
  model = NULL
  pvalue = NA
  converge = "converged"
  if (distribution == "poisson") {
    model = glm("count ~ 1", data, family = "poisson")
    if (!model$converge) {
      converge = model$converge
    }
    mu = exp(model$coefficients)
  } else if (distribution == "nb") {
    model = suppressWarnings(glm.nb("count ~ 1", data))
    size = model$theta
    if (!model$converge) {
      converge = model$converge
    }
    mu = exp(model$coefficients)
  } else if (distribution == "zinb") {
    fitResult = fitZINB(count)
    if (fitResult$converge != "converged") {
      converge = fitResult$converge
      return(cbind(pvalue, converge))
    }
    model = fitResult$modelZINB
    zeroProb = fitResult$zeroProb
    size = model$theta
    mu = exp(model$coefficients$count)
  } else {
    stop("wrong distribution specified")
  }

  if (plot) {
    if (distribution == "nb") {
      plotdist(count, distr = "nbinom", para = list(mu = mu, size = size))
    } else if (distribution == "poisson") {
      plotdist(count, distr = "pois", para = list(lambda = mu))
    } else {
      message("plot is only supported for poisson or negative binomial")
    }
  }

  # group data
  tableResult = table(count)
  value = as.numeric(names(tableResult))
  countOfValue = as.vector(tableResult)
  groupBin = combineTable(value, countOfValue, 5)

  # Check the fitting with Chi-sq Goodnees of feet test
  testStat <- c()
  testDf <- c()
  TestPvalue <- c()

  #print(i)
  observed = groupBin[3,]
  expected = numeric(length(observed))
  expected[ ] = NA
  for (j in 1:ncol(groupBin)) {
    up = groupBin[2,j]
    low = groupBin[1,j] - 1
    if (distribution == "poisson") {
      upP = ppois(up, lambda = mu)
      lowP = ppois(low, lambda = mu)
    }
    if (distribution=="nb") {
      upP = pnbinom(up, size = size, mu = mu)
      lowP = pnbinom(low, size = size, mu = mu)
    } else if (distribution == "zinb") {
      upP = pzinb(up, size = size, mu = mu, zeroProb = zeroProb)
      lowP = pzinb(low, size = size, mu = mu, zeroProb = zeroProb)
    }
    expected[j] = upP-lowP
  }
  expected[abs(expected)<1e-15] = 1e-15
  chiqTest = suppressWarnings(chisq.test(x=observed, p=expected, correct = F))
  testStat = chiqTest$statistic
  testDf = chiqTest$parameter

  testDfAdj = NULL
  if (distribution=="poisson") {
    testDfAdj <- testDf-1
  } else if (distribution=="nb") {
    testDfAdj <- testDf-2
  } else if (distribution == "zinb") {
    testDfAdj <- testDf-3
  }

  if (testDfAdj <= 0) {
    #warning("The number of bins is not enough for the test")
    converge = paste("dof is", testDfAdj)
    pvalue = NA
  } else {
    pvalue <- pchisq(testStat, df=testDfAdj, lower.tail = FALSE)
  }

  return(data.frame(pvalue, converge, stringsAsFactors = F))
}

goodnessOfFitMultiModels = function(count, dist) {
  # run goodness of fit on multiple distributions

  # input
  # count: a count vector
  # dist: a vector of distributions

  stopifnot(length(dist) >= 1)
  nModel = length(dist)
  convergeState = matrix(NA, 1, nModel)
  goodnessOfFitResult = matrix(NA, 1, nModel)
  for (j in 1:nModel) {
    result = goodnessOfFit(count, dist[j])
    convergeState[1, j] = result[, 2]
    if (result[, 2] != "converged" &&
        result[, 2] != "no zero count") {
      goodnessOfFitResult[1, j] = NA
    } else {
      goodnessOfFitResult[1, j] = result[, 1]
    }
  }

  output = cbind(goodnessOfFitResult, convergeState)
  return(output)
}

#' Run goodness of fit on a count matrix. It first downsamples all cells to a
#' common total count. Cells with total count less than this threshold will
#' be discarded
#'
#'
#' @param countMatrix a gene by cell matrix
#' @param dist a vector of distributions from "poisson", "nb", "zinb"
#' @param downsampleCount the common total count per cell used for downsampling
#' @param downsampleQuantile the quantile of the total count per cell used for
#'             downsampling. If downsampleCount is set, this will be ignored.
#'             If both are set to NULL, then no downsampling is performed
#' @param ncore the number of cores to use for parallel running
#'
#' @examples
#' \dontrun{
#' load(smallData)
#' index1 = (smallData$groupLabel == 1)
#' data1 = smallData$count[, index1]
#' group1 = smallData$groupLabel[index1]
#' goodnessResult = goodnessOfFitOnData(data1)
#' head(goodnessResult)
#' }
#'
#' @export
goodnessOfFitOnData = function(countMatrix, dist = c("nb"),
                            downsampleCount = NULL, downsampleQuantile = .1,
                            ncore = 1) {
  #browser()
  if (!is.null(downsampleCount) || !is.null(downsampleQuantile)) {
    set.seed(1234)
    downsampled = downsampleData(countMatrix,
                               groupLabel = rep(1, dim(countMatrix)[2]),
                         newTotalCounts = downsampleCount,
                         newQuantile = downsampleQuantile)
    countMatrix = downsampled[[1]]
  }

  # at least 5 cells for each gene
  kept = (rowSums(countMatrix > 0) >= 5)
  countMatrix = countMatrix[kept, , drop = F]

  nModel = length(dist)
  goodnessOfFitResult = matrix(NA, dim(countMatrix)[1], nModel)
  convergeState = matrix("", dim(countMatrix)[1], nModel)
  # browser()
  # tic = proc.time()

  # for (i in 1:dim(countMatrix)[1]) {
  #   for (j in 1:nModel) {
  #     result = goodnessOfFit(countMatrix[i, ], dist[j])
  #     convergeState[i, j] = result[, 2]
  #     if (result[, 2] != "converged" &&
  #            result[, 2] != "no zero count") {
  #         goodnessOfFitResult[i, j] = NA
  #     } else {
  #       goodnessOfFitResult[i, j] = result[, 1]
  #     }
  #   }
  #
  #   #print(result)
  # }

  # start parallel running
  library(parallel)
  maxCores = detectCores()
  if (ncore > maxCores) {
    ncore = maxCores
    print(paste("maximum available cores is", maxCores))
  }
  cluster = makeCluster(ncore)
  print(paste("number of cores used for parallel computing", ncore))

  fittedResult = parApply(cl = cluster, X = countMatrix, MARGIN = 1,
                          FUN = goodnessOfFitMultiModels, dist)
  fittedResult = t(fittedResult)

  stopCluster(cluster)
  # end parallel running

  goodnessOfFitResult = fittedResult[, 1:nModel, drop = F]
  storage.mode(goodnessOfFitResult) = "numeric"

  convergeState = fittedResult[, (nModel):(2 * nModel), drop = F]

  #browser()
  #goodnessOfFitResult[goodnessOfFitResult == -1] = 0
  geneNames = rownames(countMatrix)
  output = data.frame(geneNames, goodnessOfFitResult, convergeState,
                      stringsAsFactors = F)
  colnames(output) = c("genes", paste0("pvalue_", dist),
                       paste0("converge_", dist))

  # summarize results based on FDR
  FDR = matrix(NA, dim(countMatrix), nModel)
  for (j in 1:nModel) {
    FDR[, j] = p.adjust(goodnessOfFitResult[, j], method = "BH")
  }
  colnames(FDR) = paste0("FDR_", dist)
  output = cbind(output, FDR)

  return(output)
}
