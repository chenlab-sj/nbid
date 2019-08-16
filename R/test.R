testCovariates = function() {
  # test the handling of covariates

  set.seed(1)
  data(smallData)
  # exclude the genes with 0 in one group, because this creates large
  # differences in coefficient estimation, as well as some differences
  # in dispersions
  data = smallData$count[-c(142, 420), , drop = F]
  groups = smallData$groupLabel
  nCell = dim(data)[2]
  nCovariates = 5
  countPerCell = colSums(data)
  covariates = matrix(rnorm(nCell * nCovariates), nCell, nCovariates)
  dispMethod = "ML"

  # check without covariates
  result1 = apply(X = data, MARGIN = 1,
                     FUN = NBIDNoCov, groups,
                     countPerCell, NULL,
                     dispMethod)
  result2 = apply(X = data, MARGIN = 1,
                  FUN = NBID, groups,
                  countPerCell, NULL,
                  dispMethod)
  if (norm(t(result1)[, 1:2] - t(result2)[, 1:2], "M") > 1e-3) {
    print("FAIL")
  }

  require(MASS)
  countPerCell = colSums(data)
  LRTAndBetaNoCov = matrix(NA, dim(data)[1], 3)
  for (i in 1:dim(data)[1]) {
    resultGlmnb1 = glm.nb(data[i, ] ~ groups +
                            offset(log(countPerCell)))
    resultGlmnb0 = glm.nb(data[i, ] ~ 1 +
                            offset(log(countPerCell)))
    LRTAndBetaNoCov[i, 1] = 2 * (logLik(resultGlmnb1) - logLik(resultGlmnb0))
    LRTAndBetaNoCov[i, 2] = resultGlmnb1$coefficients[2]
    LRTAndBetaNoCov[i, 3] = 1 / resultGlmnb1$theta
  }
  # relax the threshold on differences because they are using different
  # estimations and tests
  if (norm(t(result2)[, 2, drop = F] -
           LRTAndBetaNoCov[, 2, drop = F], "M") > 0.5 ||
      norm(t(result2)[, 1, drop = F] -
           LRTAndBetaNoCov[, 1, drop = F], "M") > 2) {
    print("FAIL")
  }

  # check with covariates with glmnb
  result3 = apply(X = data, MARGIN = 1,
                  FUN = NBID, groups,
                  countPerCell, covariates,
                  dispMethod = "ML", singleDispersion = T)

  require(MASS)
  countPerCell = colSums(data)
  LRTAndBeta = matrix(NA, dim(data)[1], 3)
  for (i in 1:dim(data)[1]) {
    resultGlmnb1 = glm.nb(data[i, ] ~ groups + covariates +
                           offset(log(countPerCell)))
    resultGlmnb0 = glm.nb(data[i, ] ~ covariates +
                            offset(log(countPerCell)))
    LRTAndBeta[i, 1] = 2 * (logLik(resultGlmnb1) - logLik(resultGlmnb0))
    LRTAndBeta[i, 2] = resultGlmnb1$coefficients[2]
    LRTAndBeta[i, 3] = 1 / resultGlmnb1$theta
  }
  if (norm(t(result3)[, 2, drop = F] -
           LRTAndBeta[, 2, drop = F], "M") > 0.5 ||
      norm(t(result3)[, 1, drop = F] -
           LRTAndBeta[, 1, drop = F], "M") > 4) {
    print("FAIL")
  }

  # check with covariates with edgeR, not exactly the same, because of the
  # estimation of dispersion. edgeR uses CR, I use ML
  require(edgeR)
  groups = factor(groups)
  y = DGEList(counts = data, group = groups)
  design = model.matrix(~groups + covariates)
  y = estimateDisp(y, design, min.row.sum = 1, prior.df = 0)
  fit = glmFit(y, design, prior.count = 0)
  test = glmLRT(fit, coef=2)
  top0 = topTags(test, n = dim(data)[1])
  index = match(rownames(data), rownames(top0))
  resultEdgeR = top0@.Data[[1]][index, ]
  dispersion = y$tagwise.dispersion
  resultEdgeR = cbind(resultEdgeR, dispersion)

  print("diff in LRT between NBID and edgeR")
  diffLRT = t(result3)[, 1, drop = F] - resultEdgeR$LR
  print(summary(diffLRT))

  print("diff in log2FC between NBID and edgeR")
  diffLog2FC = t(result3)[, 2, drop = F] - resultEdgeR$logFC
  print(summary(diffLog2FC))

  index = which(abs(diffLRT) > 2)
  print(cbind(t(result3)[index, , drop = F], resultEdgeR[index, ]))

  # NBID
  result5 = apply(X = data, MARGIN = 1,
                  FUN = NBID, groups,
                  countPerCell, covariates,
                  dispMethod)

}
