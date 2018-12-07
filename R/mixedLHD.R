mixedLHD = function(
  parameters,
  indices,
  namesParameters,
  nbCondSatisfied,
  types
  ){
  
  library(lhs)
  
  integerIndices = intersect(which(types == "i"), indices)
  realIndices = intersect(which(types == "r"), indices)
  ordinalIndices =  intersect(which(types == "o"), indices)
  categoricalIndices =intersect(which(types == "c"), indices)
  
  integerNames = namesParameters[integerIndices]
  realNames <- namesParameters[realIndices]
  ordinalNames <- namesParameters[ordinalIndices]
  categoricalNames <- namesParameters[categoricalIndices]
  
  # LHD to fill
  lhd = matrix(ncol = length(indices), nrow = nbCondSatisfied)
  currentLHDColumn = 1 
  
  colNames = vector(length = ncol(lhd))
  
  if(length(categoricalNames) > 0){
    
    withCategorical = addCategorical(
                            lhd = lhd,
                            colNames = colNames, 
                            nbCondSatisfied = nbCondSatisfied,
                            categoricalNames = categoricalNames,
                            parameters = parameters,
                            currentLHDColumn = currentLHDColumn
                          )
    
    
    lhd = withCategorical[["lhd"]]
    currentLHDColumn = withCategorical[["currentLHDColumn"]]
    colNames = withCategorical[["colNames"]]
  }
  
  if(length(ordinalNames > 0)){
    
    withOrdinal = addOrdinal(
      lhd = lhd,
      colNames = colNames, 
      nbCondSatisfied = nbCondSatisfied,
      ordinalNames = ordinalNames,
      parameters = parameters,
      currentLHDColumn = currentLHDColumn
    )
    
    lhd = withOrdinal[["lhd"]]
    currentLHDColumn = withOrdinal[["currentLHDColumn"]]
    colNames = withOrdinal[["colNames"]]

  }

  if(length(integerNames) > 0){
    
    withInteger = addInteger(
      lhd = lhd,
      colNames = colNames, 
      nbCondSatisfied = nbCondSatisfied,
      integerNames = integerNames,
      parameters = parameters,
      currentLHDColumn = currentLHDColumn
    )
    
    lhd = withInteger[["lhd"]]
    currentLHDColumn = withInteger[["currentLHDColumn"]]
    colNames = withInteger[["colNames"]]
  }
  
  if(length(realNames) > 0){
    
    withReal = addReal(
      lhd = lhd,
      colNames = colNames, 
      nbCondSatisfied = nbCondSatisfied,
      realNames = realNames,
      parameters = parameters,
      currentLHDColumn = currentLHDColumn
    )
    
    lhd = withReal[["lhd"]]
    currentLHDColumn = withReal[["currentLHDColumn"]]
    colNames = withReal[["colNames"]]
    
  }
  
  return(lhd)
}

addCategorical = function(lhd, colNames, nbCondSatisfied, categoricalNames, parameters, currentLHDColumn){
  
  for (categoricalName in categoricalNames){
    
    colNames[currentLHDColumn] = categoricalName
    colnames(lhd) = colNames
    
    domain <- parameters$domain[[categoricalName]]
    
    #if the domain is larger than the number of configurations for which to sample, sample uniformly from domain
    if(length(domain) >= nbCondSatisfied){
      parameterSample = sample(domain, nbCondSatisfied, replace = FALSE)
    }
    
    #if the domain is smaller than the number of configurations for which to sample, make sure that every element is selected at least
    #floor(configurations/(length(domain)))
    #times
    else{
      
      toSample = nbCondSatisfied
      rowIndices = 1:nbCondSatisfied
      parameterSample = vector(length = nbCondSatisfied)
      
      #sample 
      while(toSample >= length(domain)){
        
        for(d in domain){
          empties = which(parameterSample == FALSE)
          
          #this distinction is necessary because sample() does not work with a vector of length 1
          if(length(empties) == 1){
            r = empties[1]
          }else{
            r = sample(empties, size = 1)
          }
          
          parameterSample[r] = d
          toSample = toSample - 1
        }
      }
      
      #sample the remainder
      if(toSample > 0){
        toInsert = sample(domain,toSample)
        
        gaps = which(parameterSample == FALSE)
        
        for(gap in gaps){
          parameterSample[gap] = toInsert[1]
          length(toInsert) = length(toInsert)-1
        }
      }
    }
    
    lhd[,currentLHDColumn] = parameterSample
    currentLHDColumn = currentLHDColumn + 1    
    
  }
  toReturn = list("lhd" = lhd, "currentLHDColumn" = currentLHDColumn, "colNames" = colNames)
  return(toReturn)
}

addOrdinal = function(lhd, colNames, nbCondSatisfied, ordinalNames, parameters, currentLHDColumn){
  
  InitialOrdinalLHD = randomLHS(n = nbCondSatisfied, k = length(ordinalNames))
  currentOrdinalLHDColumn = 1
  
  for (ordinalName in ordinalNames){
    
    colNames[currentLHDColumn] = ordinalName
    colnames(lhd) = colNames
    
    domain = parameters$domain[[ordinalName]]
    sectionSize = 1/length(domain)
    
    initialColumn = InitialOrdinalLHD[,currentOrdinalLHDColumn]
    
    parameterSample = vector()
    
    for(v in initialColumn){
      section = floor(v/sectionSize) + 1
      parameterSample = c(parameterSample, domain[section])
    }
    
    lhd[,currentLHDColumn] = parameterSample
    currentOrdinalLHDColumn = currentOrdinalLHDColumn + 1
    currentLHDColumn = currentLHDColumn + 1
  }
  
  toReturn = list("lhd" = lhd, "currentLHDColumn" = currentLHDColumn, "colNames" = colNames)
  return(toReturn)
}

addInteger = function(lhd, colNames, nbCondSatisfied, integerNames, parameters, currentLHDColumn){
  
  InitialIntegerLHD = randomLHS(n = nbCondSatisfied, k = length(integerNames))
  currentIntegerLHDColumn = 1
  
  for(integerName in integerNames){
    
    colNames[currentLHDColumn] = integerName
    colnames(lhd) = colNames
    
    parameterSample = InitialIntegerLHD[,currentIntegerLHDColumn]
    
    lhd[,currentLHDColumn] = parameterSample
    currentIntegerLHDColumn = currentIntegerLHDColumn + 1
    currentLHDColumn = currentLHDColumn + 1
  }
  
  toReturn = list("lhd" = lhd, "currentLHDColumn" = currentLHDColumn, "colNames" = colNames)
  return(toReturn)
}

addReal = function(lhd, colNames, nbCondSatisfied, realNames, parameters, currentLHDColumn){
  
  InitialRealLHD = randomLHS(n = nbCondSatisfied, k = length(realNames))
  
  currentRealLHDColumn = 1
  
  for(realName in realNames){
    
    colNames[currentLHDColumn] = realName
    colnames(lhd) = colNames
    
    parameterSample = InitialRealLHD[,currentRealLHDColumn]
    
    lhd[,currentLHDColumn] = parameterSample
    currentRealLHDColumn = currentRealLHDColumn + 1
    currentLHDColumn = currentLHDColumn + 1
  }
  
  toReturn = list("lhd" = lhd, "currentLHDColumn" = currentLHDColumn, "colNames" = colNames)
  return(toReturn)
}


fillPartialConfig = function(parameters, namesParameters, types, nbCondSatisfied,indices, configurations, digits, lhd){
  
  satisfiedCounter = 1
  
  for(idxConfiguration in seq_len(nrow(configurations))){
    configuration = configurations[idxConfiguration,]
    
    anySatisfied = FALSE
    
    for(p in indices){
      currentParameter = namesParameters[p]
      
      currentType = types[[currentParameter]]
      isSatisfied = conditionsSatisfied(parameters, configurations, currentParameter)
      anySatisfied = max(anySatisfied, isSatisfied)
      
      if(!isSatisfied){
        configuration[[p]] = NA
        next
      }
      
      else if(isFixed(currentParameter, parameters)){
        newVal = get.fixed.value(currentParameter, parameters)
      }
      
      else if(currentType == "i"){
        lowerBound = as.integer(parameters$domain[[currentParameter]][1])
        upperBound = as.integer(parameters$domain[[currentParameter]][2])
        
        normalizedVal = as.numeric(lhd[satisfiedCounter, currentParameter])
        newVal = floor(lowerBound + normalizedVal * (1 + upperBound - lowerBound))
        
        irace.assert(newVal >= lowerBound)
        irace.assert(newVal <= upperBound, eval.after = {
        sprintf("newVal = %g, normalizedVal = %g, upperBound = %g, satisfiedCounter = %g, currentParam = %s\n",
                  newVal, normalizedVal, upperBound, satisfiedCounter, currentParameter) })
        
        
      }
      
      else if(currentType == "r"){
        lowerBound <- parameters$domain[[currentParameter]][1]
        upperBound <- parameters$domain[[currentParameter]][2]
        
        normalizedVal <- as.numeric(lhd[satisfiedCounter, currentParameter])
        newVal = lowerBound + normalizedVal * (upperBound - lowerBound)
        
        #validate
        irace.assert(newVal >= lowerBound)
        irace.assert(newVal <= upperBound, eval.after = {
          sprintf("newVal = %g, normalizedVal = %g, upperBound = %g, satisfiedCounter = %g, currentParam = %s\n",
                  newVal, normalizedVal, upperBound, satisfiedCounter, currentParameter) })
      }
      
      else if(currentType %in% c("c", "o")){
        
        newVal = lhd[satisfiedCounter, currentParameter]
        
      }
      
      else{
        stop (.irace.bug.report);
      }
      
      configuration[p] = newVal
    }
    
    satisfiedCounter = satisfiedCounter + anySatisfied
    configuration = as.data.frame(configuration, stringsAsFactors = FALSE)
    configurations[idxConfiguration,] = configuration
    
  }
  
  return(configurations)
  
}