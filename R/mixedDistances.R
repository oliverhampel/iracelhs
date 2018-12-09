#returns overlap of two values (binary)    
overlap <- function(a, b){
  if(a==b){
    return(1)
  }
  else{
    return(0)
  }
  
}

#returns true if input string is "i" (integer) or "r" (real)
is_interval <- function(type){
  if (type == "i" || type == "r"){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

#returns the range for the parameter pn based on a set of configurations
get_range <- function(configurations, namesParameters, typesParameters, pn){
  
  if(is_interval(typesParameters[pn])){
    mx = max(configurations[[pn]], na.rm = TRUE)
    mn = min(configurations[[pn]], na.rm = TRUE)
    
    if (!is.null(mx) && !is.null(mn)){
      r = mx - mn
      return(r)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#returns the inverted minimal heterogeneous euclidean-overlap distance of a set of configurations
euclidean_overlap <- function(parameters, configurations){
  
  #get parameter names and types
  namesParameters <- names(parameters$conditions)
  typesParameters <- parameters$types[namesParameters]
  
  #initialize inverted min distance -> will be minimized
  inv.min.dist = NULL
  
  #calculate distance for all configurations pairwise
  for (i in  seq(1, length(configurations)-1)) {
    for (j in seq(i+1, length(configurations))){
      
      #for the curent pair of configurations calculate the distance
      d = 0
      
      for(pn in namesParameters){
        
        
        #get value for current configurations and parameter
        a = configurations[[pn]][i]
        b = configurations[[pn]][j]
        
        # add 1 to the distance if at least one of the configurations is NULL
        # calculate distance w.r.t. the current parameter and add it to the overall distance d
        if(!is.null.or.na(a) && !is.null.or.na(b)){
          if(is_interval(typesParameters[pn])){
            d = d + (abs(a-b)/get_range(configurations, namesParameters, typesParameters, pn))
          } else {
            d = d + overlap(a, b)
          }
        } else {
          #add 1 if at least one of the configurations has a NULL value for the current parameter
          d = d + 1
        }
      }
      
      
      #update inv.min.dist
      if (is.null.or.na(inv.min.dist)){
        inv.min.dist = 1/d
      } else{
        if(1/d > inv.min.dist){
          inv.min.dist = 1/d
        }
      }
    }
  }
  
  return(inv.min.dist)
  
}