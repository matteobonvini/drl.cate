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
    # print('categorical!')
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
    # print('numerical!')
    ALEPlot(X, X.model, pred.fun, j, K)
  }
}

