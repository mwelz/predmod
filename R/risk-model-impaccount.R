#' Account for imputation uncertainty in risk models
#' 
#' @param x a list of effect model objects (imputed)
#' @return A list of imputation-accounted stuff
#' 
#' @export
impaccount_risk_model <- function(x)
{
  impaccount_classcheck(x, what = "risk_model_crss")
  
  # number of imputation runs
  m <- length(x)
  
  ## prepare lists
  # benfits
  benefits <- list(absolute = NULL, relative = NULL)
  benefits$absolute <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$absolute))
  benefits$relative <- rowMeans(sapply(1:m, function(i) x[[i]]$benefits$relative))
  
  ## coefficients
  coefficients <- list(baseline = NULL, stage2 = NULL)
  
  # account for the case that baseline risk is NULL
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$coefficients$baseline ))))
  {
    brimp <- impaccount_baseline_risk(lapply(seq_len(m), function(i) x[[i]]$models$baseline))
    coefficients$baseline <- brimp$coefficients
  }
  
  # coefficients stage 2
  coefficients$stage2 <- impaccount_riskmodel_coefs(x)
  
  # risk
  risk <- list(baseline = NULL, regular = NULL, counterfactual = NULL)
  risk$regular        <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$regular)))
  risk$counterfactual <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$counterfactual)))
  
  if(all(sapply(1:m, function(i) !is.null(x[[i]]$risk$baseline ))))
  {
    risk$baseline     <- as.matrix(rowMeans(sapply(1:m, function(i) x[[i]]$risk$baseline)))
  }
  
  return(list(benefits = benefits,
              coefficients = coefficients,
              risk = risk))
} # FOR



#' internal function that accounts for imputation uncertainty the 2nd stage coefficient estimates in risk models
#' @param x a list of effect model objects (imputed)
#' @return A list of accepted and rejected coefficients
#' 
#' @noRd
impaccount_riskmodel_coefs <- function(x)
{
  ## input check not required; already done in parent function
  m <- length(x)
  
  # dimensions of the accepted models
  dim_accepted <- sapply(seq_len(m), 
                         function(i) nrow(x[[i]]$coefficients$stage2$accepted))
  
  if(identical(length(unique(dim_accepted)), 1L))
  {
    ## CASE 1: all imputation runs result in the same model
    
    decisions <- c("accepted", "rejected")
    cf_out <- list()
    
    for(decision in decisions)
    {
      p <- nrow(x[[1L]]$coefficients$stage2[[decision]])
      arr <- array(NA_real_, dim = c(p, 4L, m))
      for(i in 1:m) arr[,,i] <- x[[i]]$coefficients$stage2[[decision]]
      tmp <- impaccount_regression_array(x = arr, relative = FALSE)
      z <- tmp[, "estimate"] / tmp[, "stderr"]
      p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
      cf <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
      dimnames(cf) <- dimnames(x[[1]]$coefficients$stage2[[decision]])
      cf_out[[decision]] <- cf
    } # FOR decision
    
    # organize output
    out <- list(accepted = cf_out$accepted,
                rejected = cf_out$rejected,
                LRT = NULL)
    
  } else{
    
    ## CASE 1: some imputation runs result in different models
    
    # the full model (indicates heterogeneity) has four coefficients in
    # a cross-sectional setup, but three in a survival setup (no intercept there)
    dim_full <- ifelse(inherits(x[[1L]], what = "risk_model_crss"), 4L, 3L)
    dim_redu <- dim_full - 1L # reduced model has no interaction term
    
    
    # index in x of accepted full models 
    idx_full_acc <- which(dim_accepted == dim_full) 
    
    # index in x of rejected full models 
    idx_full_rej <- setdiff(seq_len(m), idx_full_acc) 
    
    # index in x of accepted reduced models 
    idx_redu_acc <- idx_full_rej
    
    # index in x of rejected reduced models 
    idx_redu_rej <- idx_full_acc
    
    
    ## get coefs of full and reduced models
    coefs_full_ls <- c(
      lapply(idx_full_acc, function(i) x[[i]]$coefficients$stage2$accepted),
      lapply(idx_full_rej, function(i) x[[i]]$coefficients$stage2$rejected)
    )
    
    coefs_redu_ls <- c(
      lapply(idx_redu_acc, function(i) x[[i]]$coefficients$stage2$accepted),
      lapply(idx_redu_rej, function(i) x[[i]]$coefficients$stage2$rejected)
    )
    
    
    ## make it a 3d array
    coefs_full <- array(NA_real_, dim = c(dim_full, 4L, m))
    coefs_redu <- array(NA_real_, dim = c(dim_redu, 4L, m))
    
    for(i in seq_len(m))
    {
      coefs_full[,,i] <- coefs_full_ls[[i]]
      coefs_redu[,,i] <- coefs_redu_ls[[i]]
    } # FOR m
    
    
    ## account for imputation uncertainty
    # full model
    tmp <- impaccount_regression_array(x = coefs_full, relative = FALSE)
    z <- tmp[, "estimate"] / tmp[, "stderr"]
    p <- 2.0 * stats::pnorm(abs(z), lower.tail = FALSE)
    cf_full <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
    dimnames(cf_full) <- dimnames(x[[idx_full_acc[1L]]]$coefficients$stage2$accepted)
    
    # reduced model
    tmp <- impaccount_regression_array(x = coefs_redu, relative = FALSE)
    z <- tmp[, "estimate"] / tmp[, "stderr"]
    p <- 2.0 * stats::pnorm(abs(z), lower.tail = FALSE)
    cf_redu <- cbind(tmp[, "estimate"], tmp[, "stderr"], z, p)
    dimnames(cf_redu) <- dimnames(x[[idx_redu_acc[1L]]]$coefficients$stage2$accepted)
    
    
    ## occurrence counter
    tmp <- sort(table(dim_accepted), decreasing = TRUE) # from most to least occurring
    occ_nam <- as.integer(names(tmp))
    occ_ct <- unname(tmp)
    
    winner <- ifelse(occ_nam[1L] == dim_full, "full", "reduced")
    loser <- ifelse(occ_nam[1L] == dim_full, "reduced", "full")
    
    if(winner == "full")
    {
      full_rel    <- occ_ct[1L] / m
      redu_rel    <- occ_ct[2L] / m
      cf_accepted <- cf_full
      cf_rejected <- cf_redu
    } else{
      full_rel <- occ_ct[2L] / m
      redu_rel <- occ_ct[1L] / m
      cf_accepted <- cf_redu
      cf_rejected <- cf_full
    } # IF
    
    # some information to retain
    lrt <- list(distribution = c(full = full_rel, reduced = redu_rel), 
                m = m, accepted = winner, rejected = loser)
    # print
    message(
      paste0("The LRT test has yields different models between imputation runs. Specifically, the ",
             winner, " model has been chosen in ", occ_ct[1L], " of the m=", m, " runs. ",
             "The ", loser, " model has been chosen in ", occ_ct[2L], " runs. Therefore, we treat ",
             "the ", winner, " model as the accepted model and the ", loser, " model as the rejected model.")
    )
    
    # organize output
    out <- list(accepted = cf_accepted, rejected = cf_rejected, LRT = lrt)
    
  } # IF same result across runs
  
  return(out)
} # FUN