# marginal effect
library(ALEPlot)

margin_eff <- function(X, X.model, j, pred.fun, ite, K = 40){
  # return a plot or table of marginal effect
  # X: The data frame of predictor variables to which the supervised learning model was fit.
  # X.model: The fitted supervised learning model object, typically an object to which a built-in predict command associated with that object can be applied.
  # j: index of the predictors for which the ALE plot will be calculated.
  # pred.fun: A user-supplied function that will be used to predict the response for X.model for some specified inputs
  # ite: individual treatment effects
  # K: The number of intervals into which the predictor range is divided when calculating the ALE plot effects. 
  
  if (class(X[,j]) == 'factor' | length(unique(X[,j])) <= 20){ # categorical case
    print('categorical!')
    vals <- unique(X[,j])
    bins <- length(vals)
    
    res <- data.frame(matrix(ncol = 4, nrow = bins))
    val_name <- colnames(X)[j]
    colnames(res) <- c(val_name, 'Mean', 'Upper', 'Lower')
    res[,1] <- vals

    for (i in 1:bins){
      col_val <- ite[X[,j] == vals[i]]

      # upper and lower bound
      col_val_mean <- mean(col_val)
      col_val_se <- sd(col_val)/sqrt(length(col_val))
      res[i,2] <- col_val_mean
      res[i,3] <- col_val_mean + 1.96 * col_val_se
      res[i,4] <- col_val_mean - 1.96 * col_val_se
    }
    
    res[,1] <- as.factor(res[,1])
    # a bar chart with CI
    plt <- ggplot(res, aes_string(x = val_name, y = "Mean")) +         
      geom_bar(stat="identity", color="black", width=.3, position=position_dodge()) +
      geom_errorbar(aes(ymin = Upper, ymax = Lower), width=.2,
                    position=position_dodge(.9)) +
      labs(y = 'ITEs', title = paste('Marginal effect of different values of', val_name, 'to ITEs'))
    
    print(plt)
    return(res)
    
  } else if (class(X[,j]) == 'numeric' | class(X[,j]) == 'integer'){ # continuous case
    print('numerical!')
    
    #find the vector of z values corresponding to the quantiles of X[,J]
    z= c(min(X[,j]), as.numeric(quantile(X[,j],seq(1/K,1,length.out=K), type=1)))  #vector of K+1 z values
    z = unique(z)  #necessary if X[,J] is discrete, in which case z could have repeated values 
    K = length(z)-1 #reset K to the number of unique quantile points
    fJ = numeric(K)
    #group training rows into bins based on z
    a1=as.numeric(cut(X[,j], breaks=z, include.lowest=TRUE)) #N-length index vector indicating into which z-bin the training rows fall
    X1 = X
    X2 = X
    X1[,j] = z[a1]
    X2[,j] = z[a1+1]
    y.hat1 = pred.fun(X.model=X.model, newdata = X1)
    y.hat2 = pred.fun(X.model=X.model, newdata = X2)
    Delta=y.hat2-y.hat1  #N-length vector of individual local effect values
    Delta = as.numeric(tapply(Delta, a1, mean)) #K-length vector of averaged local effect values
    fJ = c(0, cumsum(Delta)) #K+1 length vector
    #now vertically translate fJ, by subtracting its average (averaged across X[,J])
    b1 <- as.numeric(table(a1)) #frequency count of X[,J] values falling into z intervals
    fJ = fJ - sum((fJ[1:K]+fJ[2:(K+1)])/2*b1)/sum(b1)
    x <- z
    # plot(x, fJ, type="l", xlab=paste("x_",j, " (", names(X)[j], ")", sep=""), ylab= paste("f_",j,"(x_",j,")", sep=""))
    dta <- as.data.frame(cbind(x, fJ))
    ggplot(dta, aes(x = x, y = fJ)) +
      geom_line(size =0.8)+ 
      geom_point(size=1.5)+
      labs(x = eff.modif[j], y = 'ITEs', title = paste('Marginal effect of', eff.modif[j], 'to ITEs'))
  }
}

