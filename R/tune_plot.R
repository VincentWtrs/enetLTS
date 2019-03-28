tune_plot <- function(alphas, lambdas, lambda, alpha, evalCrit, index, family){
    print(paste("optimal model: lambda =", lambda, "alpha =", 
                alpha))
    lenCol <- length(alphas) * length(lambdas)
    #mycol.b <- colorRampPalette(c("black", 
    #                              "blue2", 
    #                              "purple", 
    #                              "orange", 
    #                              "yellow"))(lenCol)
    
    # NEW: REVERSING COLOR PALETTE: LIGHTER IS BETTER (Is nicer)
    mycol.b <- colorRampPalette(c("yellow", 
                                  "orange", 
                                  "purple", 
                                  "blue2", 
                                  "black"))(lenCol)
    ggmspe <- evalCrit
    rownames(ggmspe) <- lambdas
    colnames(ggmspe) <- alphas
    ggmspe <- melt(ggmspe) # NEW: Trying to force getting it from here
    if (is.null(index)){
      if (family == "binomial") {
        names(ggmspe) <- c("lambda", "alpha", "TMNLL")
        
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), y = as.factor(alpha), fill = TMNLL)) 
        + geom_tile() 
        + scale_fill_gradientn(colours = mycol.b) 
        + theme(axis.text.x = element_text(angle = -90))
        
        mspeplot <- mspeplot + ggtitle(paste0("TMNLL (minimum at lambda = ", lambda, ",alpha = ", alpha, ",  ", family, ")"))
      }
      else if (family == "gaussian") {
        names(ggmspe) <- c("lambda", "alpha", "RTMSPE")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = RTMSPE)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("RTMSPE (minimum at lambda = ", 
                                              lambda, ",alpha = ", alpha, ",  ", family, ")"))
      }
    }
    else{
      if (family == "binomial") {
        names(ggmspe) <- c("lambda", "alpha", "MNLL")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = MNLL)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("MNLL (minimum at lambda = ", 
                                              lambda, ",alpha=", alpha, ",  ", family, ")"))
      }
      else if (family == "gaussian") {
        names(ggmspe) <- c("lambda", "alpha", "RMSPE")
        mspeplot <- ggplot(ggmspe, aes(x = as.factor(lambda), 
                                       y = as.factor(alpha), fill = RMSPE)) + geom_tile() + 
          scale_fill_gradientn(colours = mycol.b) + theme(axis.text.x = element_text(angle = -90))
        mspeplot <- mspeplot + ggtitle(paste0("RMSPE (minimum at lambda=", 
                                              lambda, ",alpha = ", alpha, ",  ", family, ")"))
      }
    }
    mspeplot <- mspeplot + xlab("lambda") + ylab("alpha")
    grid.newpage() # NEW: added grid::: to force it getting the right package when changing the function in namespace # NEW NEW: REMOVED THIS
    pushViewport(viewport(layout = grid.layout(1, 1))) # From "grid" package
    print(mspeplot, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)) # From "grid" package
  }