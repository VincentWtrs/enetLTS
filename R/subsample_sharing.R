index_comparison <- function(indexall){
  ### index_comparison() FUNCTION: takes the indexall object from the warmCsteps() function and compares them to see if there is a difference
  
  # Extracting dimension sizes
  h <- dim(indexall)[1]
  length_lambda <- dim(indexall)[2]
  length_alpha <- dim(indexall)[3]
  
  # Sorting
  for(l in 1:length_lambda){
    for(a in 1:length_alpha){
      indexall[, l, a] <- sort(indexall[, l, a])
    
    }
  }
  
  index_dfrm <- vector("list", length = length_alpha)
  uniques <- vector("list", length = length_alpha)
  for(a in 1:length_alpha){
    index_dfrm[[a]] <- matrix(NA, nrow = length_lambda, ncol = h) # Rows for lambda, for each position a column
    for(l in 1:length_lambda){
      index_dfrm[[a]][l, ] <- indexall[, l, a]
    }
    # Gathering the unique best subsets
    uniques[[a]] <- unique(index_dfrm[[a]]) # Each row per list element is a unique row
  }
  
  ### ASSIGNING FOR EACH ALPHA VALUE (LIST), FOR EACH LAMBDA VALUE (POSITION WITHIN THE LIST) THE SUBSAMPLE 
  identicals <- vector("list", length = length_alpha)
  for(a in 1:length_alpha){
    identicals[[a]] <- logical(length = length_lambda)
    for(i in 1:nrow(uniques[[a]])){ # nrow not length (then it takes rowlength times rows e.g 2 x 75 = 150! (subscript out of bounds))
      for(l in 1:length_lambda){
        if(all(uniques[[a]][i, ] == index_dfrm[[a]][l, ])){
          identicals[[a]][l] <- i # We just assign the current row number
        }
      }
    }
  }
  # THIS IS IT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  return(identicals)
}
  
