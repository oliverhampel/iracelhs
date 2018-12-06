mixedLHD = function(
  parameters,
  indices,
  namesParameters,
  nbCondSatisfied,
  types = types
  ){
  
  integerIndices = intersect(which(types == "i"), indices)
  realIndices = intersect(which(types == "r"), indices)
  ordinalIndices =  intersect(which(types == "o"), indices)
  categoricalIndices =intersect(which(types == "c"), indices)
  
  integerNames = namesParameters[integerIndices]
  realNames <- namesParameters[realIndices]
  ordinalNames <- namesParameters[ordinalIndices]
  categoricalNames <- namesParameters[categoricalIndices]
  
  #generate LHD to fill
  lhd = matrix(ncol = length(indices), nrow = nbCondSatisfied)
  
  currentLHDColumn = 1 
  
  #initialize vector to store column names/parameter names
  colNames = vector(length = ncol(lhd))
  
  #iterate over categorical parameters first
  for (categoricalName in categoricalNames){
    
    #name current column of lhd
    colNames[currentLHDColumn] = categoricalName
    colnames(lhd) = colNames
    
    #get domain of current parameter
    domain <- parameters$domain[[categoricalName]]
    
    #if the domain is larger than the number of configurations for which to sample, sample uniformly from domain
    if(length(domain) >= nbCondSatisfied){

      parameterSample = sample(domain, nbCondSatisfied, replace = FALSE)
    }
    
    #if the domain is smaller than the number of configurations for which to sample, make sure that every element is selected at least
    #floor(configurations/(length(domain)))
    #times
    else{
      
      #number of values to sample
      toSample = nbCondSatisfied
      
      #create indices of the rows for which to sample
      rowIndices = 1:nbCondSatisfied
      
      #initialize vector to fill with sampled values
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
        print(gaps)
        
        for(gap in gaps){
          parameterSample[gap] = toInsert[1]
          length(toInsert) = length(toInsert)-1
        }
      }
      
    }
    
    lhd[,currentLHDColumn] = parameterSample
    currentLHDColumn = currentLHDColumn + 1    

  }  
  print(lhd)
}




