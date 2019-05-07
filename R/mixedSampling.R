mixedSampling = function(
  parameters,
  indices,
  namesParameters,
  nbCondSatisfied,
  types,
  samplingMethod
  ){

  library(lhs)
  library(randtoolbox)
  
  integerIndices = intersect(which(types == "i"), indices)
  realIndices = intersect(which(types == "r"), indices)
  ordinalIndices =  intersect(which(types == "o"), indices)
  categoricalIndices =intersect(which(types == "c"), indices)
  
  
  integerNames = namesParameters[integerIndices]
  realNames <- namesParameters[realIndices]
  ordinalNames <- namesParameters[ordinalIndices]
  categoricalNames <- namesParameters[categoricalIndices]
  
  
  #sample to fill
  sampling = matrix(ncol = length(indices), nrow = nbCondSatisfied)
  currentSamplingColumn = 1 
  
  colNames = vector(length = ncol(sampling))

  if(length(categoricalNames) > 0){
    
    withCategorical = addCategorical(
                            sampling = sampling,
                            colNames = colNames, 
                            nbCondSatisfied = nbCondSatisfied,
                            categoricalNames = categoricalNames,
                            parameters = parameters,
                            currentSamplingColumn = currentSamplingColumn
                          )
    
    
    sampling = withCategorical[["sampling"]]
    currentSamplingColumn = withCategorical[["currentSamplingColumn"]]
    colNames = withCategorical[["colNames"]]
  }
  
  nOtherParameters = length(ordinalNames) + length(integerNames) + length(realNames)
  
  if (nOtherParameters > 0) {
    
    sampling = addOthers(
      sampling = sampling,
      colNames = colNames,
      nbCondSatisfied = nbCondSatisfied,
      ordinalNames = ordinalNames,
      integerNames = integerNames,
      realNames = realNames,
      parameters = parameters,
      currentSamplingColumn = currentSamplingColumn,
      samplingMethod = samplingMethod,
      nOtherParameters = nOtherParameters
    )
  }
  
  return(sampling)
}


addCategorical = function(sampling, colNames, nbCondSatisfied, categoricalNames, parameters, currentSamplingColumn){
  
  for (categoricalName in categoricalNames){
    
    colNames[currentSamplingColumn] = categoricalName
    colnames(sampling) = colNames
    
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
      parameterSample = rep(NA, nbCondSatisfied)
      
      #sample 
      while(toSample >= length(domain)){
        
        for(d in domain){
          empties = which(is.na(parameterSample))
          
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
        
        gaps = which(is.na(parameterSample))
        
        
        library(rlist)
        
        i = 1
        for(gap in gaps){
          parameterSample[gap] = toInsert[i]
          i = i + 1
        }
      }
    }
    
    sampling[,currentSamplingColumn] = parameterSample
    currentSamplingColumn = currentSamplingColumn + 1    
    
  }
  toReturn = list("sampling" = sampling, "currentSamplingColumn" = currentSamplingColumn, "colNames" = colNames)
  return(toReturn)
}

addOthers = function(sampling, colNames, nbCondSatisfied, ordinalNames, integerNames, realNames, parameters, currentSamplingColumn, samplingMethod, nOtherParameters){
  
  if(samplingMethod == 'lhs'){
    if(nbCondSatisfied > 0){
      initialSampling = randomLHS(n = nbCondSatisfied, k = nOtherParameters)
    }else{
      initialSampling = matrix(nrow = 0, ncol = nOtherParameters)
    }
  }
  else if(samplingMethod == 'halton'){
    initialSampling = halton(n = nbCondSatisfied, dim = nOtherParameters)  
  }
  else if(samplingMethod == 'sobol'){
    initialSampling = sobol(n = nbCondSatisfied, dim = nOtherParameters)
  }
  
  currentInternalSamplingColumn = 1
  
  
  for(ordinalName in ordinalNames){
    
    colNames[currentSamplingColumn] = ordinalName
    colnames(sampling) = colNames
    
    domain = parameters$domain[[ordinalName]]
    sectionSize = 1/length(domain)
    
    if (is.null(ncol(initialSampling))){
      initialColumn = initialSampling
    }
    else if(ncol(initialSampling) >= 1){
      initialColumn = initialSampling[,currentInternalSamplingColumn] 
    }
    
    parameterSample = vector()
    
    for(v in initialColumn){
      section = floor(v/sectionSize) + 1
      parameterSample = c(parameterSample, domain[section])
    }
    
    sampling[,currentSamplingColumn] = parameterSample
    currentInternalSamplingColumn = currentInternalSamplingColumn + 1
    currentSamplingColumn = currentSamplingColumn + 1
    
  }
  
  for(integerName in integerNames){
    
    colNames[currentSamplingColumn] = integerName
    colnames(sampling) = colNames
    
    
    if (is.null(ncol(initialSampling))){
      parameterSample = initialSampling 
    }
    else if(ncol(initialSampling) >= 1){
      parameterSample = initialSampling[,currentInternalSamplingColumn]
    }
    
    
    sampling[,currentSamplingColumn] = parameterSample
    currentInternalSamplingColumn = currentInternalSamplingColumn + 1
    currentSamplingColumn = currentSamplingColumn + 1
  }
  
  for(realName in realNames){
    
    colNames[currentSamplingColumn] = realName
    colnames(sampling) = colNames
    
    if (is.null(ncol(initialSampling))){
      parameterSample = initialSampling 
    }
    else if(ncol(initialSampling) >= 1){
      parameterSample = initialSampling[,currentInternalSamplingColumn]
    }
    
    
    sampling[,currentSamplingColumn] = parameterSample
    currentInternalSamplingColumn = currentInternalSamplingColumn + 1
    currentSamplingColumn = currentSamplingColumn + 1
  }
  
  return(sampling)
  
}

fillPartialConfig = function(parameters, namesParameters, types, nbCondSatisfied, indices, configurations, digits, sampling){
  
  satisfiedCounter = 1
  
  for(idxConfiguration in seq_len(nrow(configurations))){
    configuration = configurations[idxConfiguration,]
    
    allSatisfied = FALSE
    
    for(p in indices){
      currentParameter = namesParameters[p]
      
      currentType = types[[currentParameter]]
      
      isSatisfied = conditionsSatisfied(parameters, configuration, currentParameter)
      
      allSatisfied = max(allSatisfied, isSatisfied)
      
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
      
        normalizedVal = as.numeric(sampling[satisfiedCounter, currentParameter])
        newVal = floor(lowerBound + normalizedVal * (1 + upperBound - lowerBound))
        
        irace.assert(newVal >= lowerBound)
        irace.assert(newVal <= upperBound, eval.after = {
        sprintf("newVal = %g, normalizedVal = %g, upperBound = %g, satisfiedCounter = %g, currentParam = %s\n",
                  newVal, normalizedVal, upperBound, satisfiedCounter, currentParameter) })
        
      }
      
      else if(currentType == "r"){
        lowerBound <- parameters$domain[[currentParameter]][1]
        upperBound <- parameters$domain[[currentParameter]][2]
        
        normalizedVal <- as.numeric(sampling[satisfiedCounter, currentParameter])
        newVal = lowerBound + normalizedVal * (upperBound - lowerBound)
        
        #validate
        irace.assert(newVal >= lowerBound)
        irace.assert(newVal <= upperBound, eval.after = {
          sprintf("newVal = %g, normalizedVal = %g, upperBound = %g, satisfiedCounter = %g, currentParam = %s\n",
                  newVal, normalizedVal, upperBound, satisfiedCounter, currentParameter) })
      }
      
      else if(currentType %in% c("c", "o")){
        
        newVal = sampling[satisfiedCounter, currentParameter]
        
      }
      
      else{
        stop (.irace.bug.report);
      }
      
      configuration[p] = newVal
    }
    
    satisfiedCounter = satisfiedCounter + allSatisfied
    configuration = as.data.frame(configuration, stringsAsFactors = FALSE)
    configurations[idxConfiguration,] = configuration
  }
  
  return(configurations)
  
}