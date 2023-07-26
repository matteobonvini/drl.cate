get_VIMP <- function(tau_hat, pseudo_hat, x, var.names, 
                     vimp_num_splits = 1, foldid = NULL, ...){
  # tau_hat: cate
  # pseudo_hat: pseudo outcome in stage 1
  
  lab.var.names <- unlist(lapply(var.names,
                                 function(x) gsub("[[:digit:]]", "", x[1])))
  VIM_df <- data.frame(matrix(nrow = length(lab.var.names), ncol = 3))
  rownames(VIM_df) <- unlist(lab.var.names)
  colnames(VIM_df) <- c('DR', 'Lower_Bound', 'Upper_Bound')
  
  # without sample splitting
  if (vimp_num_splits == 1){
    tau_p_hat <- mean(tau_hat)
    
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
    print('vimp sample splitting > 1')
    params <- list(...)
    
    if(is.null(foldid)) {
      s <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
    } else {
      s <- foldid
      nsplits <- length(unique(foldid))
    }
    a <- params[["a"]]
    y <- params[["y"]]
    
    for(k in vimp_num_splits){
      # in.idx <- k == s
      # ex.idx <- k != s
      
      test.idx <- k == s
      train.idx <- k != s
      
      if(all(!train.idx)) train.idx <- test.idx
      n.te <- sum(test.idx)
      n.tr <- sum(train.idx)
      
      x.tr <- x[train.idx, , drop = FALSE]
      a.tr <- a[train.idx]
      y.tr <- y[train.idx]
      x.te <- x[test.idx, , drop = FALSE]
      a.te <- a[test.idx]
      y.te <- y[test.idx]
      
      # if ((!is.matrix(v)) & (!is.data.frame(v))){
      #   v.te <- v[test.idx]
      #   v.tr <- v[train.idx]
      # } else{
      #   v.te <- v[test.idx, , drop = FALSE]
      #   v.tr <- v[train.idx, , drop = FALSE]
      # }
      
      # step 2
      pihat.vals <- option$pi.x(a = a.tr, x = x.tr, new.x = rbind(x.te, x.tr))$res
      pihat_in <- pihat.vals[1:n.te]
      pihat_ex <- pihat.vals[n.te:n.te + n.tr]
      
      mu0hat.vals <- option$mu0.x(y = y.tr, a = a.tr, x = x.tr,
                                  new.x = rbind(x.te, x.tr))$res
      mu0hat_in <- mu0hat.vals[1:n.te]
      mu0hat_ex <- mu0hat.vals[n.te:n.te + n.tr]
      
      mu1hat.vals <- option$mu1.x(y = y.tr, a = a.tr, x = x.tr,
                                  new.x = rbind(x.te, x.tr))$res
      mu1hat_in <- mu1hat.vals[1:n.te]
      mu1hat_ex <- mu1hat.vals[n.te:n.te + n.tr]
      
      pseudo_hat_in <- (a.te - pihat_in) / (pihat_in * (1 - pihat_in)) *
        (y.te - a.te * mu1hat_in - (1 - a.te) * mu0hat_in) + mu1hat_in - mu0hat_in
      pseudo_hat_ex <- (a.tr - pihat_ex) / (pihat_ex * (1 - pihat_ex)) *
        (y.tr - a.tr * mu1hat_ex - (1 - a.tr) * mu0hat_ex) + mu1hat_ex - mu0hat_ex
      
      # step3:
      drl.res <-  option$drl.x(y = pseudo, x = x.te,
                               new.x = rbind(x.te, x.tr))$res
      tau_hat_in <- drl.res[1:n.te]
      tau_hat_ex <- drl.res[n.te : n.te + n.tr]
      
      # step4:
      keep.vars <- which(!colnames(x) %in% var.names[[i]])
      fit <- drl.ite(tau_hat_in, x = x.tr[, keep.vars, drop = FALSE],
                     new.x = x.te[, keep.vars, drop = FALSE])
      tau_s_hat_in <- fit$res[, 1]
      
      # start calculating vim
      pseudo_hat_full[test.idx] <- pseudo_hat_in
      tau_hat_full[test.idx] <- tau_hat_in
      
      # need to think about how to generate multiple tau_s_hat
      # according to length(eff.modif)
      tau_s_hat_full[test.idx] <- tau_s_hat_in
    }
    tau_p_hat <- mean(tau_hat_full)
    
    # need to add for eff in eff.modif:
    theta_s_hat <- mean((pseudo_hat_full - tau_s_hat_full)^2 - (pseudo_hat_full - tau_hat_full)^2)
    theta_p_hat <- mean((pseudo_hat_full - tau_p_hat)^2 - (pseudo_hat_full - tau_hat_full)^2)
    psi_hat <- theta_s_hat / theta_p_hat
    
    phi_hat <- ((pseudo_hat_full - tau_s_hat_full)^2 -
                  psi_hat * (pseudo_hat_full - tau_p_hat)^2
                + (psi_hat - 1)*(pseudo_hat_full - tau_hat_full)^2) / theta_p_hat
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

