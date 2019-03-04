energyCriterionMixed <- function(parameters, configurations){
  namesParameters <- names(parameters$conditions)
  empty.configuration <- new.empty.configuration(parameters)
  types <- parameters$types[namesParameters]
  energyValue <- 0
  nbEnergyValues <- 0
  
  # first calculate distances between points in shared space
  sharedMask <- sapply(namesParameters,
                       isUnconditional,
                       parameters=parameters,
                       partialConfiguration=empty.configuration)
  intervalIndices <- which(types == "i" | types == "r")
  exponent <- length(intervalIndices) + 1
  sharedIntervalIndices <- intersect(intervalIndices, which(sharedMask))
  normalizedConfigs <- configurations
  
  
  
  # normalize to [-1, 1], according to Hung et al.
  for (idx in intervalIndices) {
    
    currentParam <- namesParameters[idx]
    lowerBound <- parameters$domain[[currentParam]][1]
    upperBound <- parameters$domain[[currentParam]][2]
    normalizedConfigs[,currentParam] <- normalizedConfigs[,currentParam] - 0.5 * (upperBound + lowerBound)
    normalizedConfigs[,currentParam] <- normalizedConfigs[,currentParam] / (0.5 * (upperBound - lowerBound))
  }
  
  
  # if there are shared parameters which are interval-scaled 
  # calculate energy values for them and all configurations
  if(length(sharedIntervalIndices) > 0) {
    configsSharedSubspace <- normalizedConfigs[,sharedIntervalIndices]
    if(nrow(configsSharedSubspace) > 1) {
      sharedDists <- dist(configsSharedSubspace, method="manhattan")
      energyValues <- (length(sharedIntervalIndices) / sharedDists) ** exponent
      nbEnergyValues <- nbEnergyValues + length(energyValues)
      energyValue <- energyValue + sum(energyValues)
    }
  }
  
  conditions <- parameters$conditions
  conditionStrings <- as.character(conditions)
  
  conditionStringsSet <- setdiff(conditionStrings, "TRUE")
  
  
  #iterate over all conditions
  for (conditionString in conditionStringsSet) {
    condIndices <- which(conditionStrings == conditionString)
    condIntervalIndices <- intersect(intervalIndices, condIndices)
    
    #if the current condition enables interval-scaled parameters
    if(length(condIntervalIndices) > 0) {
      #union of shared interval indices and interval indices enabled by current condition
      unionIndices <- union(sharedIntervalIndices, condIntervalIndices)
      configsSubspace <- na.omit(normalizedConfigs[,unionIndices])
      if(nrow(configsSubspace) > 1) {
        dists <- dist(configsSubspace, method="manhattan")
        energyValues <- (length(unionIndices) / dists) ** exponent
        nbEnergyValues <- nbEnergyValues + length(energyValues)
        energyValue <- energyValue + sum(energyValues)
      }
    }
  }
  energyValue <- (energyValue / nbEnergyValues) ** (1 / exponent)
  return (energyValue)
}