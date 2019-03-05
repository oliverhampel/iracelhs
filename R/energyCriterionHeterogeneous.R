energy_euclidean_overlap = function(parameters, configurations){
  
  return(energyCriterionHeterogeneous(parameters, configurations, dfunction = euclideanOverlap))
  
}


energyCriterionHeterogeneous = function(parameters, configurations, dfunction){
  
  namesParameters = names(parameters$conditions)
  empty.configuration = new.empty.configuration(parameters)
  types = parameters$types[namesParameters]
  energyValue = 0
  nbEnergyValues = 0
  
  sharedMask = sapply(namesParameters,
                       isUnconditional,
                       parameters=parameters,
                       partialConfiguration=empty.configuration)
  
  intervalIndices = which(types == "i" | types == "r")
  sharedIntervalIndices = intersect(intervalIndices, which(sharedMask))
  sharedIndices = which(sharedMask)
  
  #set criterion exponent to number of parameters + 1
  exponent = length(types) + 1
  
  #initialize normalized configs
  normalizedConfigs = configurations
  
  #we want to normalize the interval-scaled parameters to to [0;1] here because categorical and ordinal have distances ranging from 0 to 1 
  for (idx in intervalIndices){
    currentParam = namesParameters[idx]
    lowerBound = parameters$domain[[currentParam]][1]
    upperBound = parameters$domain[[currentParam]][2]
    
    normalizedConfigs[,currentParam] = (normalizedConfigs[,currentParam] - lowerBound)/(upperBound-lowerBound)
  }
  
  #calculate energy values for shared space 
  if(length(sharedIndices) > 0) {
    configsSharedSubspace = normalizedConfigs[,sharedIndices]
    if(nrow(configsSharedSubspace) > 1) {
      
      sharedDists = dfunction(configsSharedSubspace, types)
      energyValues = (length(sharedIndices)/sharedDists) ** exponent
      nbEnergyValues = nbEnergyValues + length(energyValues[!is.na(energyValues)])
      energyValue = energyValue + sum(energyValues[!is.na(energyValues)])
    }
  }
  
  conditions = parameters$conditions
  conditionStrings = as.character(conditions)
  conditionStringSet = setdiff(conditionStrings, "TRUE")
  
  
  for(conditionString in conditionStringSet){
    condIndices = which(conditionStrings == conditionString)
    
    if(length(condIndices) > 0){
      unionIndices = union(sharedIndices, condIndices)
      configsSubspace = na.omit(normalizedConfigs[,unionIndices])
      if(nrow(configsSubspace) > 1){
        dists = dfunction(configsSubspace, types)
        energyValues = (length(unionIndices)/dists) ** exponent
        nbEnergyValues = nbEnergyValues + length(energyValues[!is.na(energyValues)])
        energyValue = energyValue + sum(energyValues[!is.na(energyValues)])
      }
    }
  }
  energyValue = (energyValue/nbEnergyValues) ** (1/exponent)
  
  irace.error(energyValue)
  return(energyValue)
}