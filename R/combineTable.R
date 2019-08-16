combineTable = function(value, count, threshold) {
  # combine the table to make sure the count in each cell is above a threshold,
  # e.g., 5
  
  # input
  # value: the value 
  # count: the count of the observed value
  # threshold: a threshold value
  
  # output
  # combined: a 3 * n matrix with count >= threshold
  
  combined = NULL
  pointer = length(value)
  currentBucket = c(NA, Inf, 0) # (leftBound, rightBound, 0)
  while (pointer >= 1) {
    # add to the current bucket
    currentBucket[3] = currentBucket[3] + count[pointer]
    # check whether the bucket count >= threshold
    if (currentBucket[3] >= threshold || pointer == 1) {
      # set the left bound
      currentBucket[1] = value[pointer]
      # add to the combined
      combined = cbind(combined, currentBucket)
      # reset the currentBucket
      currentBucket = c(NA, value[pointer] - 1, 0)
    }
    pointer = pointer - 1
  }
  
  nBin = dim(combined)[2]
  
  # change the last bin to start from 0
  combined[1, nBin] = 0
  
  # make sure the last bin is >= threshold
  if (combined[3, nBin] < threshold) {
    # combine bin 1 and 2
    if (dim(combined)[2] > 1) {
      combined[1, nBin - 1] = combined[1, nBin]
      combined[3, nBin - 1] = combined[3, nBin] + combined[3, nBin - 1]
      combined = combined[, -nBin] # remove the last
    } else {
      print("not enough count with all combined")
    }
  }
  
  return(combined)
}

combineTable_test = function() {
  value = c(1, 2, 3, 5, 7, 9)
  count = c(5, 1, 5, 5, 1, 1)
  print(rbind(value, count))
  combined = combineTable(value, count, 5)
  print(combined)
}

#combineTable_test()