### Functions: Modeling Microbial Inactivation Kinetics with Disinfectant Demand
### Adapted from: Chlorine inactivation of *Elizabethkingia* spp. in water
###               https://doi.org/10.3201/eid3010.240917 
### Author: David A Holcomb
### Created: 21 June 2024
### Updated: 25 September 2024


##
### First Order Decay
##

## estimate log-linear decay rate constant k'
# argument `flat = TRUE` returns formatted estimates in dataframe, not model object
fit_loglin <- function(data,
                       formula,
                       flat = TRUE){
  
  form <- as.formula(formula)
  
  fit <- lm(form, data = data)
  
  
  if(!flat){
    
    return(fit)
    
  } else {
    
    ## fit metrics calcs  
    y <- fit$model[,1]            # original outcome data
    e_hat <- resid(fit)           # model residuals
    tss <- sum((y - mean(y))^2)   # total sum of squares
    rss <- sum((e_hat)^2)         # residual sum of squares
    r2 <- 1 - (rss / tss)         # R-squared (unreliable for NLS)
    rmse <- sqrt(mean(e_hat^2))
    dev_null <- glm(form, data = data)$null.deviance
    
    stat <- glance(fit) %>%
      mutate(tss = tss,
             dev_null = dev_null,
             r2 = r2,
             rmse = rmse) %>%
      select(sigma, rmse,
             nobs, df_resid = df.residual,
             loglik = logLik, tss,
             deviance, dev_null, dev_explain = r.squared,
             r2, aic = AIC, bic = BIC)
    
    est <- tidy(fit, conf.int = TRUE) %>%
      select(term,
             est = estimate,
             est_se = std.error,
             est_lo = conf.low,
             est_hi = conf.high,
             est_stat = statistic,
             est_pval = p.value) %>%
      bind_cols(stat)
    
    return(est)
    
  }
}


##
### Pseudo-First Order: Chick-Watson
##

## Chick-Watson rate law (no demand)
# parameters of interest (k, n) are passed log-transformed for computational reasons
law_cw <- function(c0, time, lnk, lnn){
  
  k <- exp(lnk)
  n <- exp(lnn)
  
  lnS <- -k * (c0^n) * time
  
  return(lnS)
  
}



## Chick-Watson rate law with demand
# parameters of interest (k, n) are passed log-transformed for computational reasons
law_cw_demand <- function(c0, time, kprime, lnk, lnn){
  
  k <- exp(lnk)
  n <- exp(lnn)
  
  A <- (-k * c0^n) / (n * kprime)
  B <- 1 - exp(-n * kprime * time)
  
  lnS <- A * B
  
  return(lnS)
  
}



##
### Non-Linear Log-Survival: Hom
##

## Hom model with demand
# parameters of interest (k, n, m) are passed log-transformed for computational reasons
law_hom_demand <- function(c0, time, kprime, lnk, lnn, lnm){
  
  k <- exp(lnk)
  n <- exp(lnn)
  m <- exp(lnm)
  
  A	<- (-k * m * (c0^n)) / ((n * kprime)^m)
  B <- pgamma(n * kprime * time, m)	
  
  lnS <- A * B
  
  return(lnS)
}



##
### Model Fitting & Prediction
##

### function to extract estimates and metrics from kinetic model objects
## extract parameter estimates from existing gsl_nls model object
extract_kinetic <- function(fit,
                            form_null = NULL,
                            time_var = "time"){
  
  ## R-Squared
  y <- fit$m$lhs()              # original outcome data
  e_hat <- resid(fit)           # model residuals
  tss <- sum((y - mean(y))^2)   # total sum of squares
  rss <- sum((e_hat)^2)         # residual sum of squares
  r2 <- 1 - (rss / tss)         # R-squared (unreliable for NLS)
  
  ## root mean squared error
  rmse <- sqrt(mean(e_hat^2))
  
  ## Deviance explained
  # fit null model
  if(is.null(form_null)){
    out_var <- as.character(fit$call$formula[[2]])    # outcome variable name
    form_null <- str_c(out_var, " ~ -1 + ", time_var) # nested null model formula
    message("null model not specified \nfitting linear in time with intercept fixed at 0 \n")
  }
  
  form_null <- as.formula(form_null)
  
  fit_null <- glm(form_null,
                  data = augment(fit))
  
  d_null <- deviance(fit_null)        # null deviance (residual deviance of null model)
  d_resid <- deviance(fit)            # residual deviance of proposed model
  d_explain <- 1 - (d_resid / d_null) # deviance explained
  
  
  ## standard model stats
  stat <- glance(fit) %>%
    mutate(rmse = rmse,
           tss = tss,
           r2 = r2,
           dev_explain = d_explain,
           dev_null = d_null) %>%
    select(sigma, rmse,
           nobs, df_resid = df.residual,
           loglik = logLik, tss,
           deviance, dev_null, dev_explain,
           r2, aic = AIC, bic = BIC,
           converge = isConv)
  
  est_ci <- confint(fit)
  
  ### format estimates with stats
  est <- tidy(fit) %>%
    rename(est_ln = estimate,
           est_se_ln = std.error) %>%
    mutate(est_lo_ln = est_ci[,1],
           est_hi_ln = est_ci[,2],
           term = str_remove(term, "ln"),
           est = exp(est_ln),
           est_se = exp(est_se_ln),
           est_lo = exp(est_lo_ln),
           est_hi = exp(est_hi_ln)) %>%
    bind_cols(stat) %>%
    select(term,
           est, est_se, est_lo, est_hi,
           est_stat = statistic, est_pval = p.value,
           ends_with("_ln"),
           everything())
  
  return(est)
  
}



### revised function to fit nonlinear kinetic models with nonlinear least squares
# argument `flat = TRUE` returns formatted estimates in dataframe, not model object 
fit_kinetic <- function(data,
                        formula,
                        start,
                        control = list(scale = "levenberg"),
                        flat = TRUE,
                        form_null = NULL,
                        time_var = "time"){
  
  ## ensure model formula is correctly formatted
  form <- as.formula(formula)
  
  ## fit NLS model
  fit <- gsl_nls(form,
                 data = data,
                 start = start,
                 trace = FALSE,
                 control = control)
  
  ## return model object?
  if(!flat){
    
    return(fit)
    
  } else {
    
    est <- extract_kinetic(fit,
                           form_null = form_null,
                           time_var = time_var)
    
    return(est)
    
  }
  
}


## predict log-survival from existing GSL_NLS model object
pred_kinetic <- function(fit,
                         newdata = NULL,
                         interval = "prediction",
                         level = 0.95){
  
  
  if(is.null(newdata)) newdata = augment(fit)
  
  # asymptotic predictions
  pred_raw <- predict(fit,
                      newdata = newdata,
                      interval = interval,
                      level = level) %>%
    as_tibble()
  
  pred <- bind_cols(newdata,
                    pred_raw) %>%
    rename(pred_avg = fit,
           pred_lo = lwr,
           pred_hi = upr)
  
  return(pred)
}



##
### Interpretation
##


## Calculate contact time for specified removal from model fit
# for both Hom and Chick-Watson models; set `m = 1` (default) for Chick-Watson
# assumes constant disinfectant dose
calc_time <- function(lnS, c0, k, n, m = 1){
  
  tm <- lnS / (-k * c0^n)
  
  t <- tm^(1/m)
  
  return(t)
  
}


## prepare table of CT values from model fit
tab_ct <- function(est, .dose = 0.2){
  
  cols = c(m = 1)
  
  est %>%
    select(model, nobs,
           aic, rmse, dev_explain,
           term, est) %>%
    pivot_wider(names_from = term,
                values_from = est) %>%
    tibble::add_column(., !!!cols[!names(cols) %in% names(.)]) %>%
    mutate(m = replace_na(m, 1),
           dose = .dose,
           ct_log2 = .dose * calc_time(lnS = log(0.01), c0 = .dose, k = k, n = n, m = m),
           ct_log3 = .dose * calc_time(lnS = log(0.001), c0 = .dose, k = k, n = n, m = m),
           ct_log4 = .dose * calc_time(lnS = log(0.0001), c0 = .dose, k = k, n = n, m = m),
           m = ifelse(model == "Chick-Watson", NA, m))
  
}


#
### Compare estimates
#


## Function to z-test if rate constants estimated for different data subsets are different
# rate constants and SEs should genrally be on the log-scale
test_rates <- function(k1, k1_se,
                       k2, k2_se){
  
  mu <- k1 - k2
  
  se_pool <- sqrt(k1_se^2 + k2_se^2)
  
  stat <- mu / se_pool
  
  prob <- 2 * pnorm(-abs(stat))
  
  out <- tibble(k1 = k1,
                k2 = k2,
                k_diff = mu,
                se_pool = se_pool,
                stat = stat,
                pval = prob)
  
  return(out)
  
}





#
### Full Analysis 
#

## Wrapper to conduct full analysis for given data set
analyze_kinetic <- function(data,
                            data_dis,
                            dose = 0.2){
  
  ## data
  # disinfectant decay
  df_dis <- data_dis %>%
    filter(!is.na(dis_obs)) %>%
    mutate(lnS = log( dis_obs / c0 )) %>%
    filter(time != 0)
  
  m_dis <- df_dis %>%
    fit_loglin(formula = lnS ~ -1 + time,
               flat = FALSE)
  
  est_dis <- df_dis %>%
    fit_loglin(formula = lnS ~ -1 + time) %>%
    mutate(term = ifelse(term == "time", "kprime", term),
           across(c(est, est_lo, est_hi), ~-1 * .x),
           model = "disinfectant")
  
  kprime <- est_dis %>%
    filter(term == "kprime") %>%
    select(est) %>%
    pull()
  
  # all data
  df_run <- data %>%
    mutate(kprime = kprime,
           survival = conc_nt / conc_n0,
           lnS = log(survival)) %>%
    filter(impute == "none",
           time != 0,
           !is.na(lnS))
  
  # prediction grid
  df_pred <- tibble(time = seq(0, 1, 0.01),
                    c0_dose = dose,
                    kprime = kprime)
  
  
  
  ## Fit
  # Chick-Watson
  m_cw <- fit_kinetic(df_run,
                      lnS ~ law_cw_demand(c0_dose, time, kprime, lnk, lnn),
                      start = list(lnk = 8,
                                   lnn = log(2)),
                      flat = FALSE)
  
  # Hom
  m_hom <- fit_kinetic(df_run,
                       lnS ~ law_hom_demand(c0_dose, time, kprime, lnk, lnn, lnm),
                       start = list(lnk = 20,
                                    lnn = log(2),
                                    lnm = 0.01),
                       flat = FALSE)
  
  
  ## Estimates
  est_cw <- extract_kinetic(m_cw) %>%
    mutate(model = "Chick-Watson")
  
  est_hom <- extract_kinetic(m_hom) %>%
    mutate(model = "Hom")
  
  est <- bind_rows(est_cw,
                   est_hom,
                   est_dis)
  
  
  ## Predictions
  pred_cw <- pred_kinetic(m_cw,
                          newdata = df_pred) %>%
    mutate(model = "Chick-Watson")
  
  pred_hom <- pred_kinetic(m_hom,
                           newdata = df_pred) %>%
    mutate(model = "Hom")
  
  pred_dis <- pred_kinetic(m_dis,
                           newdata = df_pred) %>%
    mutate(model = "disinfectant")
  
  pred <- bind_rows(pred_cw,
                    pred_hom,
                    pred_dis)
  
  
  
  return(list(data = df_run, est = est, pred = pred))
  
}




