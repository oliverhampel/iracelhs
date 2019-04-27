overlap = function(a, b){
  if(a == b){
    return(0)
  }else{
    return(1)
  }
}

goodall = function(a, b, clmn){
  if(a == b){
    f = sum(clmn == a)
    n = length(clmn)
    return(f*(f-1)/(n*(n-1)))
  }else{
    return(1)
  }
}

eskin2 = function(a, b, clmn){
  if(a == b){
    return(0)
  }else{
    return((4/length(unique(clmn))^2))
  }
}

occurrenceFrequency = function(a, b, clmn){
  if(a == b){
    return(0)
  }else{
    fa = sum(clmn == a)
    fb = sum(clmn == b)
    n = length(clmn)
    return(1-1/(1+(log10(n/fa)*log10(n/fb))))
  }
}

heterogeneousDistanceFunction = function(configSubspace, types, categoricalDistanceFunction){
  distmat = matrix(nrow = nrow(configSubspace), ncol = nrow(configSubspace))
  for(i in 1:(nrow(configSubspace)-1)){
    for(j in (i+1):nrow(configSubspace)){
      d = 0
      for(p in 1:ncol(configSubspace)){
        t = types[colnames(configSubspace[p])][[1]]
        if(t %in% c("c", "o" )){
          if(identical(categoricalDistanceFunction, overlap)){
            d = d + categoricalDistanceFunction(a = configSubspace[i,p], b = configSubspace[j,p])
          }else{
            d = d + categoricalDistanceFunction(a = configSubspace[i,p], b = configSubspace[j,p], clmn = configSubspace[,p])
          }
        }
        if(t %in% c("r", "i" )){
          d = d + abs(configSubspace[i,p] - configSubspace[j,p])  
        }
      }
      distmat[i, j] = d
    }
  }
  return(distmat)
}

euclideanOverlap = function(configSubspace, types){
  return(
    heterogeneousDistanceFunction(configSubspace, types, overlap)
  )
}

euclideanGoodall = function(configSubspace, types){
  return(
    heterogeneousDistanceFunction(configSubspace, types, goodall)
  )
}

euclideanEskin = function(configSubspace, types){
  return(
    heterogeneousDistanceFunction(configSubspace, types, eskin2)
  )
}

euclideanOccurrenceFrequency = function(configSubspace, types){
  return(
    heterogeneousDistanceFunction(configSubspace, types, occurrenceFrequency)
  )
}