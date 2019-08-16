#library(NBID)

checkDEAnalysis = function() {
  result = DEUsingNBID(smallData$count, smallData$groupLabel, ncore = 4)
  # check the first 10 pvalues
  # TODO: assign
  goldPvalues = c(0.902457097796305, 0.117197488526033, 0.780184466812903,
                  0.552429613512944, 0.230150047845738, 0.489924169643488,
                  0.446154932143917, 0.0315424075239331, 0.731707222282205,
                  0.152904482016356)
  if (max(abs(result[1:10, "pvalue"] - goldPvalues[1:10]) /
          goldPvalues[1:10]) > 1e-3) {
    print("FAILED on DE Analysis")
    print("gold p values")
    print(goldPvalues[1:10])
    print("calculated p values")
    print(result[1:10, "pvalue"])
  } else {
    print("PASS: DE analysis")
  }
}

checkDistribution = function() {
  index1 = (smallData$groupLabel == 1)
  data1 = smallData$count[, index1]
  group1 = smallData$groupLabel[index1]
  fitResult = checkDataDistribution(data1[1:50, ], group1, ncore = 6)
  goldAICNB = c(179.150019681101, 36.3188965084846, 30.3599533597007,
                67.0978681335691, 35.512616470689, 50.7073318750556,
                81.5034709199099, 144.506141117774, 118.946983713733,
                145.904977912094)
  if (max(abs(fitResult[1:10, "1_aicNB"] - goldAICNB[1:10]) /
          goldAICNB[1:10]) > 1e-3) {
    print("FAILED on distribution comparison")
    print("goldAICNB")
    print(goldAICNB[1:10])
    print("calculated")
    print(fitResult[1:10, "1_aicNB"])
  } else {
    print("PASS: distribution comparison")
  }

}

checkGOF = function() {
  index1 = (smallData$groupLabel == 1)
  data1 = smallData$count[, index1]
  group1 = smallData$groupLabel[index1]
  goodnessResult = goodnessOfFitOnData(data1, dist = c("poisson", "nb"),
                                       ncore = 4)

  goldValues = c(0.000170045748462803, 0.422662719867272,
                 0.112912763354193, 0.13623247841217)

  if (max(abs(goodnessResult[c(1, 4, 6, 7), 2] - goldValues) /
          goldValues) > 1e-3) {
    print("FAILED on goodness of fit")
    print("goldValues")
    print(goldAICNB)
    print("calculated")
    print(goodnessResult[c(1, 4, 6, 7), 2])
  } else {
    print("PASS: goodness of fit")
  }

}

checkMultipleGroups = function() {
  set.seed(1234)
  nGene = 30
  nGroup = 3
  nCell = 30 * nGroup
  mu = matrix(1, nGene, nCell)
  mu[11:20, 31:60] = 20
  mu[21:30, 61:90] = 20
  count = matrix(NA, nGene, nCell)
  for (i in 1:nCell) {
    count[, i] = rnbinom(n = nGene, size = 1, mu = mu[, i])
  }
  groups = c(rep(1, 30), rep(2, 30), rep(3, 30))

  index = (groups == 1 | groups == 2)
  result1vs2 = DEUsingNBID(count[, index], groups[index], ncore = 3,
                           countPerCell = 1)
  index = (groups == 1 | groups == 3)
  result1vs3 = DEUsingNBID(count[, index], groups[index], ncore = 3,
                           countPerCell = 1)
  index = (groups == 2 | groups == 3)
  result2vs3 = DEUsingNBID(count[, index], groups[index], ncore = 3,
                           countPerCell = 1)
  resultAll3 = DEUsingNBID(count, groups, ncore = 3, countPerCell = 1)
  print(cbind(result1vs2[, 1], result1vs3[, 1], result2vs3[, 1],
              resultAll3[, 1]), digits = 2)
}

main= function() {
  checkMultipleGroups()
  checkDEAnalysis()
  checkDistribution()
  checkGOF()
}

#main()
