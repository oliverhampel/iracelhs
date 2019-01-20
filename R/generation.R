#######################################
## GENERATE CONFIGURATIONS
#######################################

### FIXME: This should go within the functions that use it or in DESCRIPTION.
suppressPackageStartupMessages(library(lhs, quietly = TRUE))
suppressPackageStartupMessages(library(ParamHelpers, quietly = TRUE))

# taken from ParamHelpers:
# Convert Expressions to call (what we get from quote)
convertExpressionToCall = function(req) {
  if (is.expression(req)) {
    if (length(req) == 1) {
      return(req[[1]])
    } else {
      return(substitute(eval(x), list(x=req)))
    }
  }
  req
}

toParamSet <- function(parameters) {
  namesParameters <- names(parameters$conditions)
  params <- list()
  for (p in seq_along(namesParameters)) {
    currentParameter <- namesParameters[p]
    currentType <- parameters$types[[currentParameter]]
    condition <- parameters$conditions[[currentParameter]]
    if (convertExpressionToCall(condition) == TRUE) {
      condition <- NULL
    }
    if (currentType == "r") {
      params[[p]] <- ParamHelpers::makeNumericParam(currentParameter, 
                                      lower=parameters$domain[[currentParameter]][1],
                                      upper=parameters$domain[[currentParameter]][2],
                                      requires=condition)
    } else if (currentType == "i") {
      params[[p]] <- ParamHelpers::makeIntegerParam(currentParameter, 
                                      lower=as.integer(parameters$domain[[currentParameter]][1]),
                                      upper=as.integer(parameters$domain[[currentParameter]][2]),
                                      requires=condition)
    } else if (currentType == "c" || currentType == "o") {
      possibleValues <- parameters$domain[[currentParameter]]
      params[[p]] <- ParamHelpers::makeDiscreteParam(currentParameter, 
                                       values=possibleValues,
                                       requires=condition)
    }
  }
  ParamHelpers::makeParamSet(params=params)
}

sampleParamHelpers <- function(parameters, nbConfigurations, digits,
                           forbidden = NULL, repair = NULL){
  paramSet <- toParamSet(parameters)
  design <- ParamHelpers::generateDesign(nbConfigurations, paramSet, fun = lhs::improvedLHS)
  namesParameters <- names(parameters$conditions)
  for (p in seq_along(namesParameters)) {
    currentParameter <- namesParameters[p]
    currentType <- parameters$types[[currentParameter]]
    if (currentType == "r") {
      design[[currentParameter]] <- round(design[[currentParameter]], digits)
    }
  }
  design[[".PARENT."]] <- NA
  factorIndices <- sapply(design, is.factor)
  design[factorIndices] <- lapply(design[factorIndices], as.character)
  return(design)
}

isUnconditional <- function(parameters, partialConfiguration, paramName) {
  condition <- parameters$conditions[[paramName]]
  return(!length(all.vars(condition, max.names = 1L)))
}

## When called with an unconditional parameter, it
## must return TRUE
conditionsSatisfied <- function (parameters, partialConfiguration, paramName){
  condition <- parameters$conditions[[paramName]]
  # If there is no condition, do not waste time evaluating it.
  
  if (!length(all.vars(condition, max.names = 1L))) return(TRUE)
  
  v <- eval(condition, as.list(partialConfiguration))
  # Return TRUE if TRUE, FALSE if FALSE or NA
  v <- !is.na(v) && v
  
  return(v)
}

new.empty.configuration <- function(parameters){
  namesParameters <- names(parameters$conditions)
  newConfigurationsColnames <- c(namesParameters, ".PARENT.")
  empty.configuration <- as.list(rep(NA, length(newConfigurationsColnames)))
  names(empty.configuration) <- newConfigurationsColnames
  return(empty.configuration)
}

get.fixed.value <- function(param, parameters){
  value <- parameters$domain[[param]][1]
  type <- parameters$types[[param]]
  if (type == "i") {
    return (as.integer(value))
  } else if (type == "c" || type == "o") {
    return (value)
  } else {
    irace.assert (type == "r")
    return (as.double(value))
  }
}

energyCriterion <- function(parameters, configurations){
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
  for (conditionString in conditionStringsSet) {
    condIndices <- which(conditionStrings == conditionString)
    condIntervalIndices <- intersect(intervalIndices, condIndices)
    if(length(condIntervalIndices) > 0) {
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

correlationCriterion_old <- function (parameters, configurations){
  columnsAllNA <- apply(is.na(configurations), 2, all)
#  columnsAllNA <- as.logical(colSums(apply(configurations, 2, is.na)) == nrow(configurations))
  empty.configuration <- new.empty.configuration(parameters)
  namesParameters <- names(parameters$conditions)
  types <- parameters$types[namesParameters]

  # first calculate correlations between shared parameters
  sharedMask <- sapply(namesParameters,
                       isUnconditional,
                       parameters=parameters,
                       partialConfiguration=empty.configuration)
  indices <- which(sharedMask)
  sharedIntervalIndices <- intersect(which(types == "i" | types == "r"), indices)
  counter <- 0
  corSum <- 0
  for (i in seq_len(length(sharedIntervalIndices) - 1)) {
    for (j in (i+1):length(sharedIntervalIndices)) {
      corSum <- corSum + cor(configurations[,sharedIntervalIndices[i]],
                             configurations[,sharedIntervalIndices[j]],
                             use="complete.obs",
                             method="pearson") ** 2
      counter <- counter + 1
    }
  }

  # calculate correlations between shared and conditional parameters
  conditions <- parameters$conditions
  conditionStrings <- as.character(conditions)
  conditionStringsSet <- unique(conditionStrings)
  for (conditionString in conditionStringsSet) {
    if (conditionString != "TRUE") {
      indices <- which(conditionStrings == conditionString)
      intervalIndices <- intersect(which(types == "i" | types == "r"), indices)
      for (idx in intervalIndices) {
        if(!columnsAllNA[idx]) {
          for (sharedIdx in sharedIntervalIndices) {
            corSum <- corSum + cor(configurations[,idx],
                                   configurations[,sharedIdx],
                                   use="complete.obs",
                                   method="pearson") ** 2
            counter <- counter + 1
          }
        }
      }
    }
  }
  return (corSum / counter)
}

# alternative correlation function: calculates more pairwise correlations, but shorter code
correlationCriterion <- function (parameters, configurations){
  columnsAllNA <- apply(is.na(configurations), 2, all)
  #columnsAllNA <- as.logical(colSums(apply(configurations, 2, is.na)) == nrow(configurations))
  empty.configuration <- new.empty.configuration(parameters)
  namesParameters <- names(parameters$conditions)
  types <- parameters$types[namesParameters]

  intervalIndices <- which(types == "i" | types == "r")
  counter <- 0
  corSum <- 0
  for (i in seq_len(length(intervalIndices) - 1)) {
    for (j in (i+1):length(intervalIndices)) {
      idx1 = intervalIndices[i]
      idx2 = intervalIndices[j]
      if (!columnsAllNA[idx1] & !columnsAllNA[idx2]) {
        correlation <- cor(configurations[,idx1],
                           configurations[,idx2],
                           use="na.or.complete",
                           method="pearson") ^ 2
        if (!is.na(correlation)) {
          corSum <- corSum + correlation
          counter <- counter + 1
        }
      }
    }
  }
  return (corSum / counter)
}

weightedSum <- function(param, config) { 
  cor = correlationCriterion(param, config)
  energy = energyCriterion(param, config)
  return(energy + log10(cor))
}

both <- function(param, config) { 
  cor = correlationCriterion(param, config)
  energy = energyCriterion(param, config)
  return(c(cor, energy))
  irace.note('CURRENT PARAMS: ', param)
}

set.intervalValues <- function(lhd, nbCondSatisfied, intervalScaledIndices, intervalScaledNames){
  if(nbCondSatisfied > 1) {
    intervalValues <- (lhd - 1) / (nbCondSatisfied - 1.0)
  } else {
    # perturbed LHS (see McKay et al.) with only one point
    intervalValues <- (lhd - runif(length(intervalScaledIndices))) / nbCondSatisfied
  }
  irace.assert(all(intervalValues >= 0))
  irace.assert(all(intervalValues <= 1))
  colnames(intervalValues) <- intervalScaledNames
  return(intervalValues)
}

num.cond.satisfied <- function(parameters, indices, configurations){
  namesParameters <- names(parameters$conditions)
  satisfiedMask <- as.vector(by(configurations,
                                seq_len(nrow(configurations)),
                                conditionsSatisfied,
                                parameters = parameters,
                                paramName = namesParameters[indices[1]]))
  return(sum(satisfiedMask))
}

fillPartialConfigLegacy <- function(parameters, indices, configurations, digits,
                              intervalValues=NULL, nomialValues=NULL){
  
  #get parameter names and types
  namesParameters <- names(parameters$conditions)
  types <- parameters$types[namesParameters]
  
  #returns the number of configurations, for which the current condition is true and therefore the parameters enabled by the current condition need to be sampled 
  nbCondSatisfied <- num.cond.satisfied(parameters, indices, configurations)
  
  # if LHS for interval-scaled scaled parameters is Null (for whatever reason) generate one
  if(is.null(intervalValues)) {
    # get interval-scaled parameter indices and names of the parameters which are enabled by the current condition
    intervalScaledIndices <- intersect(which(types == "i" | types == "r"), indices)
    intervalScaledNames <- namesParameters[intervalScaledIndices]
    #generate a random LHD for all  interval-scaled parameters which are enabled by the current condition, the number of points is the number of configurations, for which the current condition is TRUE
    #values range from 1 to #points
    lhd <- randomLHD(length(intervalScaledIndices), nbCondSatisfied)
    #scale the LHD to 0..1
    intervalValues <- set.intervalValues(lhd, nbCondSatisfied, intervalScaledIndices,
                                         intervalScaledNames)
  }
  
  # if list vector for nomial-scaled parameters is Null (for whatever reason)generate one
  if(is.null(nomialValues)) {
    # get nomial parameter indices and names of the parameters which are enabled by the current condition
    nomialScaledIndices <- intersect(which(types == "c" | types == "o"), indices)
    nomialScaledNames <- namesParameters[nomialScaledIndices]
    
    #initialize a list vector, one element for each nomial scaled parameter enabled by the current condition
    nomialValues <- vector("list", length(nomialScaledNames))
    names(nomialValues) <- nomialScaledNames
    
    #determine domain for each nomial parameter enabled by the current condition and sample from it as many times as there are configurations for which the current condition is true
    for (nomialName in nomialScaledNames) {
      possibleValues <- parameters$domain[[nomialName]]
      
      #repeat value nbConSatisfied times if size of domain is 1
      if (length(possibleValues) == 1) {
        nomialValues[[nomialName]] <- rep(possibleValues, nbCondSatisfied)
      } else {
        while (length(nomialValues[[nomialName]]) < nbCondSatisfied) {
          #sample as much as possible - either size of domain or until finished
          nbSampled <- min(nbCondSatisfied - length(nomialValues[[nomialName]]),
                           length(possibleValues))
          extendedVector <- append(nomialValues[[nomialName]],
                                   sample(possibleValues, nbSampled))
          nomialValues[[nomialName]] <- extendedVector
        }
      }
    }
  }

  ## note that configurations may have more than nbCondSatisfied rows
  ## thus we need satisfiedCounter to iterate over lhs and nomialValues
  
  # nbCondSatisfied: number of configurations for which the current condition is True
  # configurations: (empty) data frame to contain the configurations
  satisfiedCounter <- 1
  #iterate over all configurations
  for (idxConfiguration in seq_len(nrow(configurations))) {
    configuration <- configurations[idxConfiguration,]
    #iterate over the indices of the parameters which are enabled/disabled by the current condition
    for (p in indices) {
      currentParameter <- namesParameters[p]
      #check if the current condition enables the current parameter for the current configuration
      isSatisfied <- conditionsSatisfied(parameters, configuration, currentParameter)
      #if the current condition is not true set to NA
      if (!isSatisfied) {
        configuration[[p]] <- NA
        next
      }
      
      # We must be careful because parameters$types does not have the
      # same order as namesParameters, because we sample in the order of the
      # conditions.
      
      #get the type of the current parameter
      currentType <- parameters$types[[currentParameter]]
      
      #if parameter is fixed - i.e. only one value is possible there is no need to sample
      if (isFixed(currentParameter, parameters)) {
        newVal <- get.fixed.value (currentParameter, parameters)
        
      # if the parameter is not fixed the value from the lhd ranged 0..1 needs to be scaled to proper parameter range
        
      # start with integer parameters  
      } 
      else if (currentType == "i") {
        #determine upper and lower bound
        lowerBound <- as.integer(parameters$domain[[currentParameter]][1])
        upperBound <- as.integer(parameters$domain[[currentParameter]][2])
        
        #scale to actual parameter range
        normalizedVal <- intervalValues[satisfiedCounter, currentParameter]
        newVal <- floor(lowerBound + normalizedVal * (1 + upperBound - lowerBound))
        # normalizedVal may be 1 so that would produce upperBound + 1:
        newVal <- min(newVal, upperBound)
        #validate newVal
        irace.assert(newVal >= lowerBound)
        irace.assert(newVal <= upperBound, eval.after = {
          sprintf("newVal = %g, normalizedVal = %g, upperBound = %g, satisfiedCounter = %g, currentParam = %s\n",
                  newVal, normalizedVal, upperBound, satisfiedCounter, currentParameter) })
        
      #real parameter
      } 
      else if (currentType == "r") {
        #determine upper and lower bound
        lowerBound <- parameters$domain[[currentParameter]][1]
        upperBound <- parameters$domain[[currentParameter]][2]

        normalizedVal <- intervalValues[satisfiedCounter, currentParameter]
        
        #validate normalized value
        irace.assert(normalizedVal >= 0)
        irace.assert(normalizedVal <= 1)
        
        #scale from 0..1 to parameter range and round 
        newVal <- lowerBound + normalizedVal * (upperBound - lowerBound)
        newVal <- round(newVal, digits)
        
        #validate
        irace.assert(newVal >= lowerBound)
        irace.assert(newVal <= upperBound, eval.after = {
          sprintf("newVal = %g, normalizedVal = %g, upperBound = %g, satisfiedCounter = %g, currentParam = %s\n",
                  newVal, normalizedVal, upperBound, satisfiedCounter, currentParameter) })
        
      #categorical and ordinal parameters
      #no need to scale, just select from nomialValues
      } 
      else if (currentType == "c" || currentType == "o") {
        newVal <- nomialValues[[currentParameter]][satisfiedCounter]
      } 
      else {
        stop (.irace.bug.report);
      }
      
      #write (scaled) value into the current configuration
      configuration[[p]] <- newVal
    }
    
    #increase satisfiedCounter if the condition is true for the current configuration - in the next iteration it needs to be sampled from the next row of the 0..1 lhs or the next element of nomialValues
    satisfiedCounter <- satisfiedCounter + isSatisfied
    
    #update/insert configurations data frame
    configuration <- as.data.frame(configuration, stringsAsFactors=FALSE)
    configurations[idxConfiguration,] <- configuration
  }
  return(configurations)
}

randomLHD <- function(dimension, numberPoints) {
  lhd <- matrix(nrow = numberPoints, ncol = dimension)
  for (i in seq_len(dimension)) {
    lhd[,i] <- sample(1:numberPoints)
  }
  return(lhd)
}

safeSample <- function(x, size, ...) {
  if(length(x) <= 1) { 
    if(!missing(size) && size == 0) {
      x[FALSE] 
    } else {
      x
    }
  } else {
    sample(x, size, ...)
  }
}

mutatedSampling <- function(sampling){
  dimension <- ncol(sampling)
  nbPoints <- nrow(sampling)
  irace.assert(dimension > 0)
  # sample number of columns to choose (1 or more) 
  nbColumns = max(1, rbinom(1, dimension, 1 / dimension))
  #sample safely nbColumn dimensions(parameters)
  chosenColumns = safeSample(seq_len(dimension), nbColumns)
  #swap two values in each selected column
  for (i in chosenColumns) {
    swapIndices <- safeSample(seq_len(nbPoints), 2)
    temp <- sampling[swapIndices[1], i]
    sampling[swapIndices[1], i] <- sampling[swapIndices[2], i]
    sampling[swapIndices[2], i] <- temp
  }
  return(sampling)
}

### Using latin hypercube sampling for the initial generation
# FIXME TODO forbidden is ignored. document/fix/remove?

sampleLHS.euclidean_overlap <- function(parameters, nbConfigurations, digits, forbidden = NULL, repair = NULL, samplingMethod){
  irace.note("Sampling ", nbConfigurations, " with LHS.euclidean_overlap\n")
  return(sampleLHS(parameters, nbConfigurations, digits, forbidden, nbEvaluations=500, objective = euclidean_overlap, legacy = FALSE, samplingMethod = samplingMethod))
}

sampleLHS.both <- function(parameters, nbConfigurations, digits, forbidden = NULL, repair = NULL, samplingMethod){

  irace.note("Sampling ", nbConfigurations, " with LHS.both\n")
  return(sampleLHS(parameters, nbConfigurations, digits, forbidden, nbEvaluations=500, objective = both, legacy = TRUE, samplingMethod = samplingMethod))
}

sampleLHS.weightedSum <- function(parameters, nbConfigurations, digits, forbidden = NULL, repair = NULL, samplingMethod){

  irace.note("Sampling ", nbConfigurations, " with LHS.weightedSum\n")

  return(sampleLHS(parameters, nbConfigurations, digits, forbidden, nbEvaluations=500, objective = weightedSum, legacy = TRUE, samplingMethod = samplingMethod))
}

sampleLHS.corr <- function(parameters, nbConfigurations, digits, forbidden = NULL, repair = NULL, samplingMethod) {
    irace.note("Sampling ", nbConfigurations, " with LHS.corr\n")
  return(sampleLHS(parameters, nbConfigurations, digits, forbidden, nbEvaluations=500, objective = correlationCriterion, legacy = TRUE, samplingMethod = samplingMethod))
}

sampleLHS.energy <- function(parameters, nbConfigurations, digits, forbidden = NULL, repair = NULL, samplingMethod) {
  irace.note("Sampling ", nbConfigurations, " with LHS.energy\n")
  return(sampleLHS(parameters, nbConfigurations, digits, forbidden, nbEvaluations=500, objective = energyCriterion, legacy = TRUE, samplingMethod = samplingMethod))
}



#samples a set of configurations using LHS, optimization criterion determined by parameter
#legacy parameter indicates whether legacy optimization criteria are used
sampleLHS <- function (parameters, nbConfigurations, digits, forbidden=NULL, nbEvaluations=1, objective=NULL, legacy, samplingMethod){
  
  #default method is energy + correlation criterion
  if(is.null(objective)) {
    objective <- function(param, config) { 
                       cor = correlationCriterion(param, config)
                       energy = energyCriterion(param, config)
                       return(c(cor, energy))
                     }
  }
  
  #read parameter names and types
  namesParameters <- names(parameters$conditions)
  types <- parameters$types[namesParameters]
  
  #initialize newConfigurations as empty data frame with NULLs, nrow = number configurations to be generated, ncol = parameters + 1
  #colnames: parameter names + .PARENT. name 
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  
  #generates an empty configuration - all parameter values and .PARENT. are NA
  empty.configuration <- new.empty.configuration(parameters)

  #get conditions as strings and a set of unique conditions
  conditions <- parameters$conditions
  conditionStrings <- as.character(conditions)
  conditionStringsSet <- unique(conditionStrings)
  
  #for each unique condition
  for (conditionString in conditionStringsSet) {
    
    #get the indices of the parameters which are enabled/disabled by the current condition
    indices <- which(conditionStrings == conditionString)
    
    #returns the number of configurations, for which the current condition is true and therefore the parameters enabled by the current condition need to be sampled 
    nbCondSatisfied <- num.cond.satisfied(parameters, indices, newConfigurations)
    
    #get interval-scaled parameter indices and names of the parameters which are enabled by the current condition
    intervalScaledIndices <- intersect(which(types == "i" | types == "r"), indices)
    intervalScaledNames <- namesParameters[intervalScaledIndices]
    
    #get nomial parameter indices and names of the parameters which are enabled by the current condition
    nomialScaledIndices <- intersect(which(types == "c" | types == "o"), indices)
    nomialScaledNames <- namesParameters[nomialScaledIndices]
    
    #if not legacy optimization generate sampling for all parameters
    if (!legacy){
      sampling = mixedSampling(
        parameters = parameters,
        indices = indices,
        namesParameters = namesParameters,
        nbCondSatisfied = nbCondSatisfied,
        types = types,
        samplingMethod = samplingMethod
      )
    }
    
    
    
    #generate a LHD for interval parameters and sample values for nomial parameters
    else{
      
      #generate a random LHD for all  interval-scaled parameters which are enabled by the current condition, the number of points is the number of configurations, for which the current condition is TRUE
      #values range from 1 to #points
      lhd <- randomLHD(length(intervalScaledIndices), nbCondSatisfied)
      
      #set column names of LHD to the names of the interval-scaled parameters which are enabled by the current condition
      colnames(lhd) <- intervalScaledNames
      
      #scale the LHD to 0..1
      intervalValues <- set.intervalValues(lhd, nbCondSatisfied, intervalScaledIndices, intervalScaledNames)
      
      #initialize a list vector, one element for each nomial scaled parameter enabled by the current condition
      nomialValues <- vector("list", length(nomialScaledNames))
      names(nomialValues) <- nomialScaledNames
      
      #determine domain for each nomial parameter enabled by the current condition and sample from it as many times as there are configurations for which the current condition is true
      for (nomialName in nomialScaledNames) {
        possibleValues <- parameters$domain[[nomialName]]
        
        #repeat value nbConSatisfied times if size of domain is 1
        if (length(possibleValues) == 1) {
          nomialValues[[nomialName]] <- rep(possibleValues, nbCondSatisfied)
        } 
        #if size of the domain is larger than 1 sample nbConSatisfied times from domain
        else {
          while (length(nomialValues[[nomialName]]) < nbCondSatisfied) {
            #sample as much as possible - either size of domain or until finished
            nbSampled <- min(nbCondSatisfied - length(nomialValues[[nomialName]]),
                             length(possibleValues))
            extendedVector <- append(nomialValues[[nomialName]],
                                     sample(possibleValues, nbSampled))
            nomialValues[[nomialName]] <- extendedVector
          }
        }
      }
    }
    
  
    if(!legacy){
      bestSampling = sampling
    }
    else{
      bestLHD = lhd
    }
    
    
    # if not legacy optimize sampling with non legacy optimization criterion 
    if(!legacy){
      
      bestCandidateConfigs = fillPartialConfig(
        parameters = parameters,
        namesParameters = namesParameters,
        types = types,
        nbCondSatisfied = nbCondSatisfied,
        indices = indices,
        configurations = newConfigurations,
        digits = digits,
        sampling = sampling
      )
      
      #optimize if there is more than one configuration which satisfies the current condition and
      #genetic algorithm has mroe than one iteration
      
      if(nbCondSatisfied > 1 && nbEvaluations > 1){
        
        bestObj = objective(parameters, bestCandidateConfigs)
        
        ### comment out because it's annoying atm
        ###print(paste(paste(indices, collapse=" "), conditionString))
        
        #optimize nbEvaluations times
        for(i in seq_len(nbEvaluations)){
          
          sampling = mutatedSampling(bestSampling)
          
          candidateConfigs = fillPartialConfig(
                                parameters = parameters, 
                                namesParameters = namesParameters, 
                                types = types,
                                nbCondSatisfied = nbCondSatisfied,
                                indices = indices, 
                                configurations = newConfigurations,
                                digits = digits,
                                sampling = sampling)
          
          obj = objective(parameters, candidateConfigs)
          
          if (all(obj <= bestObj)) {
            
            bestSampling = sampling
            bestObj = obj
            
            ### comment out because it's annoying atm
            ###irace.note('BEST: ', bestObj)
            
            bestCandidateConfigs = candidateConfigs
          
          }
          
          ### comment out because it's annoying atm
          ###cat(paste(i, obj, bestObj, sep="	", collapse="	"), "\n")
          
        }  
      }
      newConfigurations <- bestCandidateConfigs
    }
    
    #interval parameter LHD optimization
    else{
      
      #parameters of fillPartialConfigLegacy():
      # parameters:         global parameter information
      # indices:            indices of the parameters which are enabled/disabled by the current condition
      # newConfigurations:  empty data frame with NULLs, nrow = configurations, ncol = parameters + 1 (for .PARENT.)
      # digits:             precision
      # intervalValues:     scaled LHD for the interval-scaled parameters enabled by the current condition
      # nomialValues:       list-vector of configurations for the nomial-scaled parameters enabled by the current condition
      # 
      #returns a data frame containing the configurations scaled to actual parameter range and NA where condition is not True
      bestCandidateConfigs <- fillPartialConfigLegacy(parameters,
                                                      indices,
                                                      newConfigurations,
                                                      digits,
                                                      intervalValues,
                                                      nomialValues)
      
      #optimize if:
      # there is more than one interval-scaled parameter enabled by the current condition
      # and there is at least one configuration which satisfies the current condition
      # and genetic algorithm has more than one iteration
      if(length(intervalScaledIndices) > 0 && nbCondSatisfied > 1 && nbEvaluations > 1) {
        
        #evaluate objective function on initial configuration and set is as bestObj
        bestObj <- objective(parameters, bestCandidateConfigs)
        
        ### comment out because it's annoying atm
        ###print(paste(paste(indices, collapse=" "), conditionString))
        
        # optimize nbEvalutions times
        for(i in seq_len(nbEvaluations)) {
          
          #generate cofigurations scaled to actual parameter range and NA where condition is not True (just as above), but intervalValues are mutrated in each iteration
          candidateConfigs <- fillPartialConfigLegacy(parameters,
                                                      indices,
                                                      newConfigurations,
                                                      digits,
                                                      intervalValues,
                                                      nomialValues)
          
          #evaluate objective function  on current configurations 
          obj <- objective(parameters, candidateConfigs)
          
          #if all objectives are smaller or equal currently best objective update best objective and update best LHD
          #this code minimizes the objective
          
          if (all(obj <= bestObj)) {
            bestObj <- obj
            
            ### comment out because it's annoying atm
            ###irace.note('BEST: ', bestObj)
            
            #best LHD ist kept for mutation
            bestLHD <- lhd
            #best configurations are stored
            bestCandidateConfigs <- candidateConfigs
          }
          
          ### comment out because it's annoying atm
          ###cat(paste(i, obj, bestObj, sep="	", collapse="	"), "\n")
          
          # mutate for next iteration
          lhd <- mutatedSampling(bestLHD)
          
          #scale the LHD to 0..1 and verify to get new intervalValues
          intervalValues <- (lhd - 1) / (nbCondSatisfied - 1.0)
          irace.assert(all(intervalValues >= 0))
          irace.assert(all(intervalValues <= 1))
        }
      }
      
      newConfigurations <- bestCandidateConfigs
    }
  }
  return (newConfigurations)
}


### Uniform sampling for the initial generation
sampleUniform <- function (parameters, nbConfigurations, digits,
                           forbidden = NULL, repair = NULL){
  # Sample new configurations.
  if (getOption(".irace.debug.level", default = 0) >= 1) {
    irace.note("Sampling ", nbConfigurations,
               " configurations from uniform distribution\n")
  }
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  empty.configuration <- new.empty.configuration(parameters)

  for (idxConfiguration in seq_len(nbConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      configuration <- empty.configuration
      for (p in seq_along(namesParameters)) {
        currentParameter <- namesParameters[p]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          configuration[[p]] <- NA
          next
        }
        # FIXME: We must be careful because parameters$types does not have the
        # same order as namesParameters, because we sample in the order of the
        # conditions.
        currentType <- parameters$types[[currentParameter]]
        if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType == "i") {
          lowerBound <- as.integer(parameters$domain[[currentParameter]][1])
          upperBound <- as.integer(parameters$domain[[currentParameter]][2])
          newVal <- floor(runif(1, min = lowerBound, max = 1 + upperBound))
        } else if (currentType == "r") {
          lowerBound <- parameters$domain[[currentParameter]][1]
          upperBound <- parameters$domain[[currentParameter]][2]
          newVal <- runif(1, as.double(lowerBound), as.double(upperBound))
          newVal <- round(newVal, digits)
        } else if (currentType == "c" || currentType == "o") {
          possibleValues <- parameters$domain[[currentParameter]]
          newVal <- sample(possibleValues, 1)
        } else {
          stop (.irace.bug.report);
        }
        configuration[[p]] <- newVal
      }
      configuration <- as.data.frame(configuration, stringsAsFactors=FALSE)
      if (!is.null(repair)) {
        configuration <- repair(configuration, parameters, digits)
      }

      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}

# To be called the first time before the second race (with indexIter =
# 2) Nb configurations is the number of configurations at the end
# included the elite ones obtained from the previous iteration
sampleModel <- function (parameters, eliteConfigurations, model,
                         nbNewConfigurations, digits, forbidden = NULL,
                         repair = NULL){
  if (nbNewConfigurations <= 0) {
    irace.error ("The number of configurations to generate appears to be negative or zero.")
  }
  namesParameters <- names(parameters$conditions)
  newConfigurations  <-
    as.data.frame(matrix(nrow = nbNewConfigurations,
                         ncol = length(namesParameters) + 1,
                         dimnames = list(NULL, c(namesParameters, ".PARENT."))
                         ))
  empty.configuration <- new.empty.configuration(parameters)

  for (idxConfiguration in seq_len(nbNewConfigurations)) {
    forbidden.retries <- 0
    while (forbidden.retries < 100) {
      # Choose the elite which will be the parent.
      indexEliteParent <- sample.int (n = nrow(eliteConfigurations), size = 1,
                                      prob = eliteConfigurations[[".WEIGHT."]])
      eliteParent <- eliteConfigurations[indexEliteParent, ]
      idEliteParent <- eliteParent[[".ID."]]
      configuration <- empty.configuration
      configuration[[".PARENT."]] <- idEliteParent

      # Sample a value for every parameter of the new configuration.
      for (p in seq_along(namesParameters)) {
        # FIXME: We must be careful because parameters$types does not
        # have the same order as parameters$conditions. Ideally, we
        # should fix this or make it impossible to confuse them.
        currentParameter <- namesParameters[p]
        currentType <- parameters$types[[currentParameter]]
        if (!conditionsSatisfied(parameters, configuration, currentParameter)) {
          # Some conditions are unsatisfied.
          # Should be useless, NA is ?always? assigned when matrix created
          newVal <- NA

        } else if (isFixed(currentParameter, parameters)) {
          # We don't even need to sample, there is only one possible value !
          newVal <- get.fixed.value (currentParameter, parameters)
          # The parameter is not a fixed and should be sampled
        } else if (currentType == "i" || currentType == "r") {
          lowerBound <- paramLowerBound(currentParameter, parameters)
          upperBound <- paramUpperBound(currentParameter, parameters)
          mean <- as.numeric(eliteParent[currentParameter])
          if (is.na(mean)) {
            # The elite parent does not have any value for this
            # parameter, let's sample uniformly.
            newVal <- ifelse(currentType == "i",
                             floor(runif(1, min = lowerBound, max = 1 + upperBound)),
                             runif(1, lowerBound, upperBound))
          } else {
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]]
            newVal <- ifelse(currentType == "i",
                             rtnorm(1, mean + 0.5, stdDev, lowerBound, upperBound + 1) - 0.5,
                             rtnorm(1, mean, stdDev, lowerBound, upperBound))
          }
          newVal <- ifelse(currentType == "i", round(newVal),
                           round(newVal, digits))

        } else if (currentType == "o") {
          possibleValues <- paramDomain(currentParameter, parameters)
          value <- eliteParent[currentParameter]

          if (is.na(value)) {
            # The elite parent does not have any value for this
            # parameter, let's sample uniformly
            newVal <- sample(possibleValues, 1)
          } else {
            # Find the position within the vector of possible
            # values to determine the equivalent integer.
            mean <- match(value, possibleValues) # Return index of value in array
            stdDev <- model[[currentParameter]][[as.character(idEliteParent)]]

            # Sample with truncated normal distribution as an integer.
            newValAsInt <- round(rtnorm(1, mean + 0.5, stdDev, 1,
                                        length(possibleValues) + 1) - 0.5)

            # Get back to categorical values, find the one corresponding to the
            # newVal
            newVal <- possibleValues[newValAsInt]
          }
        } else if (currentType == "c") {
          # FIXME: Why is idEliteParent character?
          # FIXME: Why the model is <parameter><Parent>? It makes more sense to be <Parent><parameter>.
          probVector <- model[[currentParameter]][[as.character(idEliteParent)]]
          possibleValues <- paramDomain(currentParameter, parameters)
          newVal <- sample(x = possibleValues, size = 1, prob = probVector)
        } else {
          stop (.irace.bug.report)
        }
        configuration[[p]] <- newVal
      }

      configuration <- as.data.frame(configuration, stringsAsFactors=FALSE)
      if (!is.null(repair)) {
        configuration <- repair(configuration, parameters, digits)
      }
      if (is.null(forbidden)
          || nrow(checkForbidden(configuration, forbidden)) == 1) {
        newConfigurations[idxConfiguration,] <- configuration
        break
      }
      forbidden.retries <- forbidden.retries + 1
    }
    if (forbidden.retries >= 100) {
      irace.error("irace tried 100 times to sample from the model a configuration not forbidden without success, perhaps your constraints are too strict?")
    }
  }
  return (newConfigurations)
}
