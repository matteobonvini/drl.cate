get_VIMP <- function(tau_hat, pseudo_hat, x, var.names, num_split = 1){
  # tau_hat: cate
  # pseudo_hat: pseudo outcome in stage 1
  
  # without sample splitting
  if (num_split == 1){
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
  }
  
  # with sample splitting
  else {
    # todo
  }
  return(VIM_df)
}

