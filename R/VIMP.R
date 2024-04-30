#' get_vimp
#'
#' This function calculates Variable IMPortance measures following the procedure
#' in Hines et al (2022) https://arxiv.org/pdf/2204.06030.pdf
#'
#' @param cate.fit object from fitting the cate using drl.cate
#'
#' @export
get_vimp <- function(cate.fit, var.names, lab.var.names,
                     x = NULL, y = NULL, a = NULL,
                     nsplits = NULL, pi.x = NULL, mu1.x = NULL, mu0.x = NULL,
                     drl.x = NULL){
  VIM_df <- data.frame(matrix(nrow=length(var.names), ncol=5))
  rownames(VIM_df) <- lab.var.names
  colnames(VIM_df) <- c('psi', "theta_s", "theta_p", 'lb', 'ub')

  if(!is.null(cate.fit)){
    x <- cate.fit$x
    n <- nrow(x)
    drl.x <- cate.fit$drl.x
    tau_hat <- cate.fit$cate.x.res[["cate.x.sample"]][["dr"]][, 1]
    pseudo_hat <- cate.fit$cate.x.res[["pseudo"]][["dr"]]
    ites.x.tr <- cate.fit$cate.x.res[["cate.x.sample.tr"]][["dr"]]
    foldid <- cate.fit$foldid
    nsplits <- length(unique(foldid))
  } else {
    n <- length(y)
    foldid <- sample(rep(1:nsplits, ceiling(n / nsplits))[1:n])
    pseudo_hat <- tau_hat <- rep(NA, n)
    ites.x.tr <- vector("list", nsplits)

    for(k in 1:nsplits){
      # in.idx <- k == s
      # ex.idx <- k != s

      test.idx <- k == foldid
      train.idx <- k != foldid

      if(all(!train.idx)) train.idx <- test.idx
      n.te <- sum(test.idx)
      n.tr <- sum(train.idx)

      x.tr <- x[train.idx, , drop = FALSE]
      a.tr <- a[train.idx]
      y.tr <- y[train.idx]
      x.te <- x[test.idx, , drop = FALSE]
      a.te <- a[test.idx]
      y.te <- y[test.idx]
      print('initializing x y a successfully!')
      print(a)
      # step 2
      pihat.vals <- pi.x(a = a.tr, x = x.tr, new.x = rbind(x.te, x.tr))$res

      print('calculating pi successfully!')

      pihat_in <- pihat.vals[1:n.te]
      pihat_ex <- pihat.vals[-c(1:n.te)]

      mu0hat.vals <- mu0.x(y = y.tr, a = a.tr, x = x.tr,
                                  new.x = rbind(x.te, x.tr))$res
      mu0hat_in <- mu0hat.vals[1:n.te]
      mu0hat_ex <- mu0hat.vals[-c(1:n.te)]

      mu1hat.vals <- mu1.x(y = y.tr, a = a.tr, x = x.tr,
                           new.x = rbind(x.te, x.tr))$res
      mu1hat_in <- mu1hat.vals[1:n.te]
      mu1hat_ex <- mu1hat.vals[-c(1:n.te)]
      print('calculating muhat successfully!')

      # Step 3 calculate pseudo outcomes and tau_hat
      pseudo_hat_in <- (a.te - pihat_in) / (pihat_in * (1 - pihat_in)) *
        (y.te - a.te * mu1hat_in - (1 - a.te) * mu0hat_in) + mu1hat_in - mu0hat_in
      pseudo_hat_ex <- (a.tr - pihat_ex) / (pihat_ex * (1 - pihat_ex)) *
        (y.tr - a.tr * mu1hat_ex - (1 - a.tr) * mu0hat_ex) + mu1hat_ex - mu0hat_ex

      pseudo_hat[test.idx] <- pseudo_hat_in

      drl.vals.x <- drl.x(pseudo = pseudo_hat_ex, x = x.tr,
                          new.x = rbind(x.te, x.tr))$res
      tau_hat[test.idx] <- drl.vals.x[1:n.te, 1]
      ites.x.tr[[k]] <- drl.vals.x[-c(1:n.te), 1]

    }
  }

    tau_p_hat <- mean(tau_hat)
    theta_p_hat <- mean((pseudo_hat - tau_p_hat)^2) - mean((pseudo_hat - tau_hat)^2)

    for(i in 1:length(var.names)){
      keep.vars <- which(!colnames(x) %in% var.names[[i]])
      tau_s_hat <- rep(NA, n)
      for(k in 1:nsplits){
        idx.te <- k == foldid
        if(length(unique(foldid)) > 1) idx.tr <- k != foldid
        else idx.tr <- idx.te

        fit <- drl.x(pseudo = ites.x.tr[[k]],
                     x = x[idx.tr, keep.vars, drop = FALSE],
                     new.x = x[idx.te, keep.vars, drop = FALSE])
        tau_s_hat[idx.te] <- fit$res[, 1]

      }

      theta_s_hat <- mean((pseudo_hat - tau_s_hat)^2) - mean((pseudo_hat - tau_hat)^2)

      psi_hat <- theta_s_hat / theta_p_hat

      phi_hat <- ((pseudo_hat - tau_s_hat)^2 -
                    psi_hat * (pseudo_hat - tau_p_hat)^2
                  + (psi_hat - 1)*(pseudo_hat - tau_hat)^2) / theta_p_hat
      var.val <- 1/length(phi_hat)^2 * sum((phi_hat)^2)

      VIM_df[i, "psi"] <- psi_hat
      VIM_df[i, "theta_s"] <- theta_s_hat
      VIM_df[i, "theta_p"] <- theta_p_hat
      VIM_df[i, "lb"] <- max(VIM_df[i, 1] - 1.96 * sqrt(var.val), 0)
      VIM_df[i, "ub"] <- min(VIM_df[i, 1] + 1.96 * sqrt(var.val), 1)
      # print(VIM_df[1:i, ])
    }
  return(VIM_df)
}


#' draw_VIMP
#'
#'Function to draw the VIMP results
#'
#' @param VIMP_2b a get_vimp object
#'
#' @export
draw_VIMP <- function(VIM_2b){
  VIM_2b <- cbind(as.matrix(row.names(VIM_2b), nrow(VIM_2b), 1), VIM_2b)
  colnames(VIM_2b) <- c('variable', 'psi', 'lb', 'ub')
  VIM_2b <- VIM_2b[order(VIM_2b$DR, decreasing = FALSE),]
  order_2b <- VIM_2b$variable

  fig_2b <- ggplot(VIM_2b, aes(x = factor(variable, level = order_2b), y = DR))+
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Lower_Bound, ymax = Upper_Bound), width=.1, position=position_dodge(0.2))+
    theme(axis.title.y=element_blank(),
          axis.text.y = element_text(color = "grey20", size = 10),
          axis.text.x = element_text(color = "grey20", size = 10)) +
    coord_flip()
  fig_2b
}

