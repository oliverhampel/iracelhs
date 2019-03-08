euclideanOverlap = function(configSubspace, types){
  distmat = matrix(nrow = nrow(configSubspace), ncol = nrow(configSubspace))
  for(i in 1:(nrow(configSubspace)-1)){
    for(j in (i+1):nrow(configSubspace)){
      d = 0
      for(p in 1:ncol(configSubspace)){
        t = types[colnames(configSubspace[p])][[1]]
        if(t %in% c("c", "o" )){
            d = d + overlap(configSubspace[i,p], configSubspace[j,p])
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

overlap = function(a, b){
  if(a == b){
    return(0)
  }else{
    return(1)
  }
}

euclideanGoodall = function(configSubspace, types){
  distmat = matrix(nrow = nrow(configSubspace), ncol = nrow(configSubspace))
  for(i in 1:(nrow(configSubspace)-1)){
    for(j in (i+1):nrow(configSubspace)){
    d = 0
      for(p in 1:ncol(configSubspace)){
        t = types[colnames(configSubspace[p])][[1]]
        if(t %in% c("c", "o" )){
          d = d + goodall(configSubspace[i,p], configSubspace[j,p], clmn = configSubspace[,p])
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

goodall = function(a, b, clmn){
  if(a == b){
    f = sum(clmn == a)
    n = length(clmn)
    return(f*(f-1)/(n*(n-1)))
  }else{
    return(1)
  }
}

euclideanEskin = function(configSubspace, types){
  distmat = matrix(nrow = nrow(configSubspace), ncol = nrow(configSubspace))
  for(i in 1:(nrow(configSubspace)-1)){
    for(j in (i+1):nrow(configSubspace)){
      d = 0
      for(p in 1:ncol(configSubspace)){
        t = types[colnames(configSubspace[p])][[1]]
        if(t %in% c("c", "o" )){
          d = d + eskin2(configSubspace[i,p], configSubspace[j,p], clmn = configSubspace[,p])
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

eskin2 = function(a, b, clmn){
  if(a == b){
    return(0)
  }else{
    return((4/length(unique(clmn))^2))
  }
}

euclideanOccurrenceFrequency = function(configSubspace, types){
  distmat = matrix(nrow = nrow(configSubspace), ncol = nrow(configSubspace))
  for(i in 1:(nrow(configSubspace)-1)){
    for(j in (i+1):nrow(configSubspace)){
      d = 0
      for(p in 1:ncol(configSubspace)){
        t = types[colnames(configSubspace[p])][[1]]
        if(t %in% c("c", "o" )){
          d = d + occurrenceFrequency(configSubspace[i,p], configSubspace[j,p], clmn = configSubspace[,p])
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