# library(dplyr)
# library(survey)
# library(tableone)
# library(glmnet)
# source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/DrawTableOneFxn.R')
# source('~/Desktop/EHR-mining/Scripts/AnalysisScripts/ProcessTableOneFxn.R')

# df_og <- df

library(dplyr)
library(survival)

getDataSet <- function(filter_var) {
   if (!filter_var %in% c('iDAP', 'allDAP', 'i1DAP', 'i2DAP', 'sDAP')) {
      cat('Filtering variable not recognized.')
      break
   }
   
   df <- read.table(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/FULL_MRSA_BACT.txt',  # change to working directory
                    header = TRUE,
                    sep = '\t')
   df <- tibble(df)
   
   
   ##### START FILTERING ##### 
   # iDAP i1DAP i2DAP sDAP allDAP
   if (filter_var == 'sDAP') {
      df <- df[df$TRT %in% c('VAN', 'sDAP'), ]
      
   } else if (filter_var == 'i2DAP') {
      df <- df[df$TRT %in% c('VAN', 'iDAP') | (df$TRT == 'sDAP' & df$DAP_SWITCH_TIME < 2), ]
      
   } else if (filter_var == 'i1DAP') {
      df <- df[df$TRT %in% c('VAN', 'iDAP') | (df$TRT == 'sDAP' & df$DAP_SWITCH_TIME < 1), ]
      
   } else if (filter_var == 'iDAP') {
      df <- df[df$TRT %in% c('VAN', 'iDAP'), ]
   }
   df$TRT <- as.integer(df$TRT == 'VAN')
   df <- df[names(df) != 'DAP_SWITCH_TIME']
   ##### END FILTERING #####
   
   
   ##### START DEAL WITH MISSING LAB VALUES #####
   # leave impute presence
   w_cols <- grep('^PRE_', names(df), value=TRUE)
   df[w_cols] <- sapply(df[w_cols],
                        FUN = function(vec) {
                           w <- which(is.na(vec))
                           mean_val <- mean(vec, na.rm=T)
                           vec[w] <- mean_val
                           return(vec)
                        })
   
   df <- df %>%
      mutate(CATEGORY = case_when(
         CATEGORY == 1 ~ 'academic',
         CATEGORY == 2 ~ 'rural',
         CATEGORY == 3 ~ 'community',
         CATEGORY == 4 ~ 'regional'
      )) %>%
      mutate(academic = as.integer(CATEGORY == 'academic'),
             community = as.integer(CATEGORY == 'community'),
             regional = as.integer(CATEGORY == 'regional')) %>%
      select(-CATEGORY)
   bin_vars <- names(df)[sapply(df, function(x) length(unique(x)) == 2L)]
   df <- df %>%
      mutate(across(.cols = !!bin_vars,
                    .fns = ~ as.factor(.)),
             status = ifelse(time == 31, 0, 1))
   
   return(df)
}

getFormula <- function(set) {
   base_formula <- 'Surv(time, status) ~ TRT'
   if (set == 'univariate') {                         return(base_formula)
   } else if (set == 'multivariate (all)') {          vars <- names(df)[!names(df) %in% c('TRT', 'time', 'status')]
   } else if (set == 'multivariate (trt - 0.01)') {   vars <- trt_vars_01
   } else if (set == 'multivariate (trt - 0.05)') {   vars <- trt_vars_05
   } else if (set == 'multivariate (trt - 0.1)') {    vars <- trt_vars_1
   } else if (set == 'multivariate (time - 0.01)') {  vars <- surv_vars_01
   } else if (set == 'multivariate (time - 0.05)') {  vars <- surv_vars_05
   } else if (set == 'multivariate (time - 0.1)') {   vars <- surv_vars_1
   } else if (set == 'multivariate (both - 0.01)') {  vars <- intersect(trt_vars_01, surv_vars_01)
   } else if (set == 'multivariate (both - 0.05)') {  vars <- intersect(trt_vars_05, surv_vars_05)
   } else if (set == 'multivariate (both - 0.1)') {   vars <- intersect(trt_vars_1, surv_vars_1) }
   if (length(vars) == 0L) return(NULL)
   formula <- paste0(base_formula, ' + ', paste(vars, collapse=' + '))
   return(formula)
}

parseSurvModel <- function(s) {
   s <- summary(s)$coefficients['TRT1', ]
   pval <- unname(s['Pr(>|z|)'])
   est <- unname(s['exp(coef)'])
   se <- unname(s['se(coef)'])
   return(c('pval' = pval, 'est' = est, 'se' = se))
}



filter_vars <- c('iDAP', 'allDAP', 'i1DAP', 'i2DAP', 'sDAP')
resDF <- setNames(vector(mode = 'list', 
                         length = length(filter_vars)),
                  filter_vars)
filter_var <- 'sDAP'

for (filter_var in filter_vars) {
   df <- getDataSet(filter_var)
   cat(filter_var, nrow(df), '\n')
   
   
   # define groups of covariates
   all_vars <- names(df)[!names(df) %in% c('TRT', 'time', 'status')]
   
   ps_mod <- glm(as.formula(paste0('TRT ~ ', paste(names(df)[!names(df) %in% c('TRT', 'time', 'status')], collapse=' + '))),
                 family = binomial(),
                 data = df)
   sig_vars <- summary(ps_mod)$coefficients
   sig_vars <- sig_vars[rownames(sig_vars) != '(Intercept)', ]
   rownames(sig_vars) <- gsub('1$', '', rownames(sig_vars))
   trt_vars_01 <- rownames(sig_vars)[sig_vars[, 'Pr(>|z|)'] < 0.01]
   trt_vars_05 <- rownames(sig_vars)[sig_vars[, 'Pr(>|z|)'] < 0.05]
   trt_vars_1 <- rownames(sig_vars)[sig_vars[, 'Pr(>|z|)'] < 0.1]
   
   surv_mod <- coxph(as.formula(paste0('Surv(time, status) ~ ', paste(names(df)[!names(df) %in% c('TRT', 'time', 'status')], collapse=' + '))),
                     data = df) %>% summary()
   surv_vars <- surv_mod$coefficients
   surv_vars[is.na(surv_vars)] <- 1
   rownames(surv_vars) <- gsub('1$', '', rownames(surv_vars))
   surv_vars_01 <- rownames(surv_vars)[surv_vars[, 'Pr(>|z|)'] < 0.01]
   surv_vars_05 <- rownames(surv_vars)[surv_vars[, 'Pr(>|z|)'] < 0.05]
   surv_vars_1 <- rownames(surv_vars)[surv_vars[, 'Pr(>|z|)'] < 0.1]
   
   
   # define all relevant adjustment sets
   adj_set <- c('univariate',
                'multivariate (all)', 
                'multivariate (trt - 0.01)', 'multivariate (trt - 0.05)', 'multivariate (trt - 0.1)',
                'multivariate (time - 0.01)', 'multivariate (time - 0.05)', 'multivariate (time - 0.1)', 
                'multivariate (both - 0.01)', 'multivariate (both - 0.05)', 'multivariate (both - 0.1)')
   
   
   # Let's add propensity scores and IPT weights: 
   marg_prob_trt <- sum(df$TRT == 1) / nrow(df)
   ps_df <- tibble(
      score = predict(ps_mod, type = 'response'),
      TRT = ps_mod$model$TRT
   ) %>%
      mutate(trimmed_score = case_when(
         score < 0.01 ~ 0.01,
         score > 0.99 ~ 0.99,
         .default = score
      )) %>%
      mutate(weights = case_when(
         TRT == 1L ~ 1 / trimmed_score,
         TRT == 0L ~ 1 / (1 - trimmed_score)
      )) %>%
      # A * (1 - p) + (1 - A) * p
      mutate(overlap_weights = case_when(
         TRT == 1L ~ 1 - score,
         TRT == 0L ~ score
      )) %>%
      mutate(stab_weights = case_when(
         TRT == 1L ~ marg_prob_trt / score,
         TRT == 0L ~ (1 - marg_prob_trt) / (1 - score)
      ))
   cat('Pseudo population sizes:\n')
   cat('IPTW:', sum(ps_df$weights) / nrow(df), '\n')
   cat('OW:', sum(ps_df$overlap_weights) / nrow(df), '\n')
   cat('SW:', sum(ps_df$stab_weights) / nrow(df), '\n')
   cat('Fraction extreme propensity scores:', sum(ps_df$score < 0.01 | ps_df$score > 0.99) / nrow(ps_df), '\n')
   
   # par(mfrow=c(2,1))
   # plot(NA, xlim=c(0, 1), ylim=c(0,6), ylab='Density', xlab='Propensity score')
   # lines(density(ps_df$score[ps_df$TRT == 1], bw='SJ'), col='blue', lwd=2)
   # lines(density(ps_df$score[ps_df$TRT == 0], bw='SJ'), col='red', lwd=2)
   # 
   # #par(mfrow=c(1,1))
   # plot(NA, xlim=c(0, 1), ylim=c(0,5.5), ylab='Density', xlab='Propensity score')
   # lines(density(ps_df$overlap_weights[ps_df$TRT == 1], bw='SJ'), col='blue', lwd=2)
   # lines(density(ps_df$overlap_weights[ps_df$TRT == 0], bw='SJ'), col='red', lwd=2)
   # 
   # #par(mfrow=c(1,1))
   # plot(NA, xlim=c(1,34), ylim=c(0,5.5), ylab='Density', xlab='Propensity score')
   # lines(density(ps_df$weights[ps_df$TRT == 1], bw='SJ'), col='blue', lwd=2)
   # lines(density(ps_df$weights[ps_df$TRT == 0], bw='SJ'), col='red', lwd=2)
   
   X <- df[, !names(df) %in% c('TRT', 'time', 'status')]
   binvars <- X %>% select(where(is.factor)) %>% names
   X <- as.matrix(X %>% mutate(across(where(is.factor), ~ as.integer(.) - 1)))
   A <- as.integer(df$TRT)-1
   bal <- col_w_smd(mat=X, treat=A, weights=ps_df$stab_weights, bin.vars=binvars) %>% round(5)
   if (any(bal != 0)) {
      cat('Balance failed!\n')
      #break
   }
   
   # lower_ps_limit <- max(tapply(ps_df$score, ps_df$TRT, min))
   # upper_ps_limit <- min(tapply(ps_df$score, ps_df$TRT, max))
   # ps_df %>% filter(score <= upper_ps_limit, score >= lower_ps_limit) %>% count(TRT)
   # ps_df %>% count(TRT)
   # cat('', sum(ps_df$score >= lower_ps_limit & ps_df$score <= upper_ps_limit) / nrow(df), '\n')
   
   
   # run regressions
   results <- data.frame(
      pval = rep(NA_real_, length(adj_set) * 4),
      est = rep(NA_real_, length(adj_set) * 4),
      se = rep(NA_real_, length(adj_set) * 4),
      desc = paste0(rep(adj_set, each=4), c('', ' - IPTW', ' - OW', ' - SW'))
   )
   for (i in seq_along(adj_set)) {
      set <- adj_set[i]
      formula <- getFormula(set)
      if (is.null(formula)) {
         # cat(i, '')
         next
      }
      # unweighted
      res <- df %>% coxph(formula = as.formula(formula), data = .) %>% parseSurvModel()
      results[i * 4 - 3, 1:3] <- res
      # IPTW - ATE
      res <- df %>% coxph(formula = as.formula(formula), data = ., weights = ps_df$weights) %>% parseSurvModel()
      results[i * 4 - 2, 1:3] <- res
      # OW - ATO
      res <- df %>% coxph(formula = as.formula(formula), data = ., weights = ps_df$overlap_weights) %>% parseSurvModel()
      results[i * 4 - 1, 1:3] <- res
      # stabilized weights - AT??
      res <- df %>% coxph(formula = as.formula(formula), data = ., weights = ps_df$stab_weights) %>% parseSurvModel()
      results[i * 4, 1:3] <- res
   }
   # results <- results[grep('([^W]|SW)$', results$desc), ]
   resDF[[filter_var]] <- results
   
   w <- grep('W$', results$desc, invert=TRUE)
   ypos <- seq_along(w)
   yshift <- 0.2
   putArrows <- function(w, col='black', yshift=0) {
      points(x=results$est[w], y=ypos+yshift, pch=16, col=col)
      arrows(x0=results$est[w] - qnorm(0.975) * results$se[w],
             x1=results$est[w] + qnorm(0.975) * results$se[w],
             y0=seq_along(w) + yshift,
             angle=90, length=0.05, code=3, col=col)
   }

   par(mfrow=c(1,1), mar=c(4,11,2,1), mgp=c(1.5, 0.5, 0), tck=-0.015)
   plot(NA, xlim=c(0.5, 2.5), ylim=range(ypos) + c(-0.5, 0.5), ylab='', xlab='Hazard ratio', yaxt='n', log='x', xaxt='n', main=filter_var)
   axis(side=1, at=c(0.5, 1, 1.5, 2))
   abline(v = 1, lty=2)
   #points(x=results$est[w], y=ypos, pch=16)
   #points(x=results$est[-w], y=ypos-yshift, pch=16, col='red')
   putArrows(w=w)
   n_extra_rows <- nrow(results) / length(adj_set)
   for (i in 2:n_extra_rows)
      putArrows(w=seq(i, nrow(results), n_extra_rows), col='red', yshift=(i-1)*yshift)
   axis(side=2, at=seq_along(w), labels=results$desc[w], las=1)
   
   
   # remove unnecessary variables
   rm(df, all_vars, ps_mod, sig_vars, trt_vars_01, trt_vars_05, trt_vars_1, surv_mod, surv_vars, surv_vars_01, surv_vars_05, surv_vars_1,
      set, formula, i, res, results, adj_set, ps_df)
   cat('\n')
}















# Create table 1
tOg <- CreateTableOne(vars=names(df)[!names(df) %in% c('TRT', 'time')], strata='TRT', data=df)
drawTableOne(processTableOne(tOg))







# matching using MatchIt
# optimal matching - becuase nearest didn't balance
psm_mod <- MatchIt::matchit(formula = as.formula(paste0('TRT ~ ', paste(grep('TRT|time', names(df), invert=TRUE, value=TRUE), collapse=' + '))),
                   method = 'optimal',
                   data = df)
df_match <- MatchIt::match.data(psm_mod)

tOg <- tableone::CreateTableOne(vars=names(df_match)[!names(df_match) %in% c('TRT', 'time', 'distance', 'weights', 'subclass')], strata='TRT', data=df_match, smd=TRUE)
drawTableOne(processTableOne(tOg))

df_match <- df_match %>% mutate(status = ifelse(time == 31, 0, 1))
df_match %>% coxph(Surv(time, status) ~ TRT, data=.) %>% summary() %>% getTrtCoef() # 0.146
df_match %>%
   coxph(formula = as.formula(paste0('Surv(time, status) ~ TRT + ',
                                     paste(names(df_match)[!names(df_match) %in% c('TRT', 'time', 'distance', 'weights', 'subclass', 'status')],
                                           collapse=' + '))), 
         data=.) %>% summary() %>% getTrtCoef() # 0.016
df_match <- df_match %>% select(-status)



# inverse PS weighting
# use propensity scores in ps_df
#df$weights <- 
ps_df <- ps_df %>%
   mutate(weights = case_when(
      TRT == 1L ~ 1 / score,
      TRT == 0L ~ 1 / (1 - score)
   ))
# df$weights <- ps_df$weights

# weight
par(mfrow=c(1,1))
plot(NA, xlim=c(0, 34), ylim=c(0,5.5), ylab='Density', xlab='Weights')
lines(density(ps_df$weights[ps_df$TRT == 1], bw='SJ'), col='blue', lwd=2)
lines(density(ps_df$weights[ps_df$TRT == 0], bw='SJ'), col='red', lwd=2)

# original
tOg <- tableone::CreateTableOne(vars=names(df)[!names(df) %in% c('TRT', 'time', 'status')], strata='TRT', data=df)
drawTableOne(processTableOne(tOg))

# weighted
survey_design <- svydesign(ids = ~ 1, data=df, weights = ~ ps_df$weights)
tOg_wtd <- svyCreateTableOne(vars=names(df)[!names(df) %in% c('TRT', 'time', 'status')], strata='TRT', data=survey_design)
drawTableOne(processTableOne(tOg_wtd), add=TRUE)

df <- df %>% mutate(status = ifelse(time == 31, 0, 1))
df %>% coxph(Surv(time, status) ~ TRT, data=.) %>% summary() # 0.0076
df %>%
   coxph(formula = as.formula(paste0('Surv(time, status) ~ TRT + ',
                                     paste(names(df)[!names(df) %in% c('TRT', 'time', 'subclass', 'status')],
                                           collapse=' + '))), 
         data=.) %>% summary() # 0.0015
df %>%
   coxph(formula = as.formula(paste0('Surv(time, status) ~ TRT + ',
                                     paste(names(df)[!names(df) %in% c('TRT', 'time', 'subclass', 'status')],
                                           collapse=' + '))), 
         data=., 
         weights=ps_df$weights) %>% summary() # 0.0065
sig_vars <- rownames(summary(ps_mod)$coefficients)[summary(ps_mod)$coefficients[, 'Pr(>|z|)'] < 0.05]
sig_vars <- sig_vars[sig_vars != '(Intercept)']
sig_vars <- gsub('1$', '', sig_vars)
df %>%
   coxph(formula = as.formula(paste0('Surv(time, status) ~ TRT + ', paste(sig_vars, collapse=' + '))),
         data = .,
         weights = ps_df$weights) %>% summary() # 0.012
sig_vars <- rownames(summary(ps_mod)$coefficients)[summary(ps_mod)$coefficients[, 'Pr(>|z|)'] < 0.1]
sig_vars <- sig_vars[sig_vars != '(Intercept)']
sig_vars <- gsub('1$', '', sig_vars)
df %>%
   coxph(formula = as.formula(paste0('Surv(time, status) ~ TRT + ', paste(sig_vars, collapse=' + '))),
         data = .,
         weights = ps_df$weights) %>% summary() # 0.027
df <- df %>% select(-status)






## Propensity score calculation with regularization??
library(glmnet)

fit <- glmnet(x = df %>% select(!c(TRT, time)), y = df$TRT, family = 'binomial')
plot(fit, label=TRUE)
print(fit)
coef(fit, s=0.01)



# 







as.formula(paste0('Surv(time, status) ~ TRT + ', paste(vars, collapse=' + ')))
# setup plot area
ORplot <- function(trt) {
   par(mar=c(5, 4, 2.75, 3.5))
   plot(NA, xlim=c(0.5,4.5), ylim=c(1/5, 5), log='y', xaxt='n', xlab='', ylab='', yaxt='n', main='30-day cox proportional hazards', cex.main=cex)
   title(ylab='Hazard ratio', line=2.5, cex.lab=cex)
   title(xlab='Model', line=3.5, cex.lab=cex)
   axis(side=2, at=c(1/4, 1/2, 1, 2, 3, 4), labels=as.character(c(1/4, 1/2, 1, 2, 3, 4)), las=1)
   text(x=1:2, y=1/6, labels=c('Unadjusted', 'Inverse propensity\nscore weighted'), adj=c(0.5, 1), xpd=NA, cex=cex)
   abline(h = 1, lty=2)
   arrows(x0=0.5, y0=1/4.1, y1=4.1, length=0.1, code=3)
   text(x=0.45, y=c(1/4.8, 4.8), adj=0, labels=paste0('favors ', c(trt[1], trt[2])), cex=cex)
}
ORplot(trt)

s <- summary(coxph(formula = formula, data=df))
est <- s$coefficients[1, 'coef']
se <- qnorm(0.975) * s$coefficients[1, 'se(coef)']
points(x=1, y=exp(est), pch=16, cex=cex)
arrows(x0=1, y0=exp(est-se), y1=exp(est+se), code=3, angle=90, length=0.05, lwd=cex)

prop_model <- df %>% # predicts trt1
   glm(formula = paste0('trt1 ~ ', paste(vars, collapse=' + ')),
       family = binomial(),
       data = .)

df$prop_score <- predict.glm(prop_model, newdata=df, type='response')
df <- df %>%
   mutate(prop_weights = case_when(
      trt1 == 1L ~ 1 / prop_score,
      trt1 == 0L ~ 1 / (1 - prop_score)
   ))

s <- summary(coxph(formula = formula, data=df, weights=prop_weights))
est <- s$coefficients['trt1', 'coef']
se <- qnorm(0.975) * s$coefficients[1, 'se(coef)']
points(x=2, y=exp(est), pch=16, cex=cex)
arrows(x0=2, y0=exp(est-se), y1=exp(est+se), code=3, angle=90, length=0.05, lwd=cex)

coefs_lin_raw <- glm(formula=paste0('time ~ trt1 + ', paste(vars, collapse=' + ')), data=df) %>% summary()
coefs_lin_adj <- glm(formula=paste0('time ~ trt1 + ', paste(vars, collapse=' + ')), data=df, weights=prop_weights) %>% summary()
print(coefs_lin_raw$coefficients['trt1', 'Pr(>|t|)'])
print(coefs_lin_adj$coefficients['trt1', 'Pr(>|t|)'])

coefs_log_raw <- df %>% mutate(d30 = time < 30) %>% glm(formula=paste0('d30 ~ trt1 + ', paste(vars, collapse=' + ')), data=., family=binomial()) %>% summary()
coefs_log_adj <- df %>% mutate(d30 = time < 30) %>% glm(formula=paste0('d30 ~ trt1 + ', paste(vars, collapse=' + ')), data=., family=binomial(), weights=prop_weights) %>% summary()
print(coefs_log_raw$coefficients['trt1', 'Pr(>|z|)'])
print(coefs_log_adj$coefficients['trt1', 'Pr(>|z|)'])

df <- df %>% select(-prop_score, -prop_weights, -trt1)



