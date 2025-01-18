library(survival)
library(dplyr)
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/MRSA_BSI_DAP_VAN.Rdata')
load(file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/LabValues_DF_MRSA.Rdata')

df <- dfx; rm(dfx)
vars <- names(df)[!names(df) %in% c('ORDER_DATE', 'FIRST_VAN_TIME', 'FIRST_DAP_TIME',
                                    'FACILITY', 'FACILITY_TRANSFER', 'TRTs', 'DAYS_SINCE_PRV', 'VAN_MIC_2', 'VAN_MIC',
                                    grep('^(POST|PRE)_', names(df), value=TRUE, ignore.case=TRUE), 
                                    'CATEGORY_TRANSFER', 'RESULT_DAY', 'RESULT_DATE', 'RECENT_MRSA',
                                    'FIRST_ABX_DAY_REL_ORDER', 'FIRST_ABX_TIME_REL_ORDER', 'LENGTH_OF_STAY', 'TRANSFER_TIME',
                                    'RECUR', 'ABX_START_TIME', 'ABX_END_TIME', 'ABX_START_DAY', 'ABX_END_DAY', 'TRANSFER', 'year')]
df <- df %>%
   select(!!vars) %>%
   left_join(
      x = .,
      y = lab_vals,
      by = join_by(PERSON_ID, ORDER_DAY)
   ) %>%
   select(-PERSON_ID, -ORDER_DAY, -time) %>%
   rename(time = time_censored) %>%
   mutate(
      NumOtherIsolates = as.numeric(NumOtherIsolates),
      across(.cols = where(is.logical),
             .fns = ~ as.integer(.)),
      CATEGORY = as.numeric(case_when(
         CATEGORY == 'academic' ~ 1,
         CATEGORY == 'rural' ~ 2,
         CATEGORY == 'community' ~ 3,
         CATEGORY == 'regional' ~ 4
      ))
   )

write.table(x = df, 
            file = '~/Desktop/EHR/EHR work/RdataFiles/causal_prep/MRSA_bacteremia/FULL_MRSA_BACT.txt',
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)


##### MULTIVARIATE REGRESSION ####
# Survival
s <- df %>%
   filter(TRT != 'sDAP') %>%
   mutate(across(contains('PRE'), ~ tidyr::replace_na(., mean(., na.rm=T)))) %>%
   select(!c(TRT, time, PERSON_ID, ORDER_DAY)) %>%
   coxph(Surv(time_censored, status) ~ ., data=.) %>%
   summary()
s <- s$coefficients
vars_surv <- rownames(s)[s[, 'Pr(>|z|)'] < 0.1]

# Treatment
t <- df %>%
   filter(TRT != 'sDAP') %>%
   mutate(across(contains('PRE'), ~ tidyr::replace_na(., mean(., na.rm=T)))) %>%
   select(!c(time, time_censored, status, PERSON_ID, ORDER_DAY)) %>%
   mutate(TRT = as.factor(TRT)) %>%
   glm(TRT ~ ., family=binomial(link='logit'), data=.) %>%
   summary()
t <- t$coefficients
vars_trt <- rownames(t)[t[, 'Pr(>|z|)'] < 0.1]

# Overalp
length(vars_surv) # 21
length(vars_trt) # 26
length(intersect(vars_surv, vars_trt)) # 4 LOL
##### END #####



select_vars <- c('or1', 'or5', 'everything', 'and1', 'and5')
filter_vars <- c('iDAP', 'allDAP', 'i1DAP', 'i2DAP', 'sDAP')
handle_missing_vars <- c('impute', 'presence')#, 'leave')
include_outcome_vars <- c('yes', 'no')
search_alg_vars <- c('pc', 'fges', 'grasp-fci', 'fci', 'boss', 'grasp', 'gfci', 'rfci')
search_alpha_vars <- c(0.05, 0.01, 0.1)
params <- expand.grid(filter_var=filter_vars, 
                      select_var=select_vars, 
                      handle_missing_var=handle_missing_vars, 
                      include_outcome_var=include_outcome_vars, 
                      search_alg_var=search_alg_vars, 
                      search_alpha_var=search_alpha_vars)

# for score-based algorithms, they all have every possible p-value, but that is unecessary
# remove the p-value and then take unique
score_based_algs <- c('boss', 'fges')
params$search_alpha_var[params$search_alg_var %in% score_based_algs] <- NA
params <- unique(params)
rownames(params) <- seq_len(nrow(params))
nrow(params)

# programatically construct the .dag file
write.table(x = params[1:3,], 
            file = '~/Desktop/CausalDiscovery/job_files/search.map', 
            quote = FALSE,
            col.names = FALSE,
            row.names = TRUE,
            sep = ' ')






