#to write down some comments

draw_VIMP <- function(vimp_df){
  # a helper function to illustrate vimp
  
  vimp_df <- cbind(as.matrix(row.names(vimp_df), nrow(vimp_df), 1), vimp_df)
  colnames(vimp_df) <- c('variable', 'DR', 'Lower_Bound', 'Upper_Bound')
  vimp_df <- vimp_df[order(vimp_df$DR, decreasing = FALSE),]
  order_2b <- vimp_df$variable
  
  fig_2b <- ggplot(vimp_df, aes(x = factor(variable, level = order_2b), y = DR))+
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width=.1, position=position_dodge(0.2))+
    theme(axis.title.y=element_blank(), 
          axis.text.y = element_text(color = "grey20", size = 10),
          axis.text.x = element_text(color = "grey20", size = 10)) +
    ylim(-0.1, 1) +
    coord_flip() 
  fig_2b
}

# the  main function
get_VIMP <- function(tau_hat, pseudo_hat, x, var.names){
  #
  #
  #
  # var.names: variable of interest, i.e. effect modifier
  tau_p_hat <- mean(tau_hat)
  lab.var.names <- unlist(lapply(var.names, 
                                 function(x) gsub("[[:digit:]]", "", x[1])))
  
  VIM_df <- data.frame(matrix(nrow = length(lab.var.names), ncol = 3))
  rownames(VIM_df) <- unlist(lab.var.names)
  colnames(VIM_df) <- c('DR', 'Lower_Bound', 'Upper_Bound')
  # x.mat <- data[, xnames, drop = FALSE]
  # v.mat <- data[, eff.modif, drop = FALSE]
  
  for(i in 1:length(var.names)){
    keep.vars <- which(!colnames(x) %in% var.names[[i]])
    fit <- drl.ite(tau_hat, x = x[, keep.vars, drop = FALSE],
                   new.x = x[, keep.vars, drop = FALSE])
    tau_s_hat <- fit$res[, 1]
    theta_s_hat <- mean((pseudo_hat - tau_s_hat)^2 - (pseudo_hat - tau_hat)^2)
    theta_p_hat <- mean((pseudo_hat - tau_p_hat)^2 - (pseudo_hat - tau_hat)^2)
    
    psi_hat <- theta_s_hat / theta_p_hat
    
    phi_hat <- ((pseudo_hat - tau_s_hat)^2 - 
                  psi_hat * (pseudo_hat - tau_p_hat)^2 
                + (psi_hat - 1)*(pseudo_hat - tau_hat)^2) / theta_p_hat
    var.val <- 1/length(phi_hat)^2 * sum((phi_hat)^2)
    
    VIM_df[i, 1] <- psi_hat
    # VIM_3b[i, 2] <- VIM_3b[i, 1] - 1.96 * sqrt(var)
    # VIM_3b[i, 3] <- VIM_3b[i, 1] + 1.96 * sqrt(var)
    VIM_df[i, 2] <- max(VIM_df[i, 1] - 1.96 * sqrt(var.val), 0)
    VIM_df[i, 3] <- min(VIM_df[i, 1] + 1.96 * sqrt(var.val), 1)
    
    print(c(VIM_df[i,1], VIM_df[i,2], VIM_df[i,3],theta_s_hat, theta_p_hat))
  }
  return(VIM_df)
}

