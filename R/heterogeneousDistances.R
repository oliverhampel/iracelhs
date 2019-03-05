euclideanOverlap = function(configSubspace, types){
  
  
  distmat = matrix(nrow = nrow(configSubspace), ncol = nrow(configSubspace))
  
  for(i in 1:(nrow(configSubspace)-1)){
    
    for(j in (i+1):nrow(configSubspace)){
      
      d = 0
      
      for(p in 1:ncol(configSubspace)){
        
        t = types[colnames(configSubspace[p])][[1]]
        
        if(t %in% c("c", "o" )){
          if(configSubspace[i,p]!=configSubspace[j,p]){
            d = d + 1
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