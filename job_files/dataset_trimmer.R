# This script takes as input the full MRSA bacteremia EHR dataset
# and several arguments deciding how to filter and select variables
# the only external package required is `survival`

# get args from command line
args <- commandArgs(trailingOnly = TRUE)
cat('COMMAND LINE FLAGS AND ARGS:\n', args, '\n\n')
flags <- gsub('--', '', args[seq(1, length(args), 2)])
vars <- args[seq(2, length(args), 2)]

# assign command line arguments to variables
job_num <- vars[flags == 'job_num']
filter_var <- vars[flags == 'filter_var']
select_var <- vars[flags == 'select_var']
handle_missing_var <- vars[flags == 'handle_missing_var']
include_outcome_var <- vars[flags == 'include_outcome_var']

# check the output to ensure everything is kosher
cat('job_num', job_num, typeof(job_num), '\n')
cat('filter_var', filter_var, typeof(filter_var), '\n')
cat('select_var', select_var, typeof(select_var), '\n')
cat('handle_missing_var', handle_missing_var, typeof(handle_missing_var), '\n')
cat('include_outcome_var', include_outcome_var, typeof(include_outcome_var), '\n\n')

# job_num <- 1
# filter_var <- 'iDAP'
# select_var <- 'and1'
# handle_missing_var <- 'impute'
# include_outcome_var <- 'yes'

# read in FULL dataset from home directory
# this is the base dataset that will be passed to each worker node
# next, filter cases and select variables according to input params:
#     filter_vars
#     handle_missing
#     select_vars
#     include_outcome
library(survival)
df <- read.table(file = 'FULL_MRSA_BACT.txt',  # change to working directory
                 header = TRUE,
                 sep = '\t')


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
if (handle_missing_var == 'impute') {
   df[w_cols] <- sapply(df[w_cols],
                        FUN = function(vec) {
                           w <- which(is.na(vec))
                           mean_val <- mean(vec, na.rm=T)
                           vec[w] <- mean_val
                        })
   
} else if (handle_missing_var == 'presence') {
   df[w_cols] <- sapply(df[w_cols], function(vec) as.integer(!is.na(vec)))
   # if we do this, some variables (labs) will be PERFECTLY or NEARLY PERFECTLY correlated
   # we need to account for this by removing all but 1 or combining them
   marked_for_removal <- character()
   for (v1 in seq_along(w_cols)[-length(w_cols)]) {
      vec1 <- df[w_cols[v1]]
      for (v2 in seq_along(w_cols)[(v1+1):length(w_cols)]) {
         vec2 <- df[w_cols[v2]]
         if (cor(vec1, vec2) == 1L) {
            marked_for_removal <- c(marked_for_removal, w_cols[v2])
         }
      }
   }
   marked_for_removal <- unique(marked_for_removal)
   df <- df[!names(df) %in% marked_for_removal]
}
##### END DEAL WITH MISSING LAB VALUES #####


##### START VARIABLE SELECTION #####
covars <- names(df)[!names(df) %in% c('TRT', 'time', 'status')]
if (select_var == 'everything') { # only if we are not keeping everything
   keep_vars <- covars
   
} else {
   varsDF <- data.frame(
      var = covars, 
      pval_surv = numeric(length(covars)), 
      pval_trt = numeric(length(covars))
   )
   for (v in seq_along(covars)) {
      # Survival
      s <- summary(coxph(formula = as.formula(paste0('Surv(time, status) ~ ', covars[v])), data=df))
      varsDF$pval_surv[v] <- s$coefficients[,'Pr(>|z|)']
      
      # Treatment
      t <- summary(glm(formula = as.formula(paste0('TRT ~ ', covars[v])), family=binomial(link = 'logit'), data=df))
      t <- t(as.matrix(t$coefficients[-1,]))
      varsDF$pval_trt[v] <- t[, 'Pr(>|z|)']
   }
   
   # everything or1 or5 and1 and5
   alpha <- unname(c('or1' = 0.1, 'or5' = 0.05, 'and1' = 0.1, 'and5' = 0.05)[select_var])
   if (grepl('or', select_var)) {
      w <- which(varsDF$pval_surv < alpha | varsDF$pval_trt < alpha)
      
   } else {
      w <- which(varsDF$pval_surv < alpha & varsDF$pval_trt < alpha)
   }
   keep_vars <- varsDF$var[w]
}

# only keep necessary variables for causal search
df <- df[c('TRT', 'time', keep_vars)]
##### END VARIABLE SELECTION #####


# include outcome?
if (include_outcome_var == 'no') {
   df <- df[names(df) != 'time']
}


##### CREATE KNOWLEDGE FILE #####
sink('knowledge.txt')
cat('/knowledge\n\naddtemporal\n')

temporal_count <- 1
if (any(keep_vars == 'AGE')) {
   cat(temporal_count, 'AGE\n')
   temporal_count <- temporal_count + 1
}

cat(temporal_count, ' ', paste(keep_vars, collapse=' '), '\n', sep='')
temporal_count <- temporal_count + 1

cat(temporal_count, 'TRT\n')
temporal_count <- temporal_count + 1

if (include_outcome_var == 'yes') {
   cat(temporal_count, 'time\n')
   temporal_count <- temporal_count + 1
}
cat('\nforbiddirect\n')
cat('\nrequiredirect\n')
sink()
##### END #####


##### DUMP LIST OF VARIABLES TO FILE #####
write.table(x = names(df), 
            file = paste0('variables_', job_num, '.txt'),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
##### END #####



# print out values to ensure correctness
cat(nrow(df), ncol(df), filter_var, select_var, handle_missing_var, include_outcome_var, '\n\nR SCRIPT DONE RUNNING\nWRITE TRIMMED DATA TO TXT FILE AND RUN CAUSAL SEARCH\n\n')


# write trimmed data set to (worker node?) directory
write.table(x = df,
            file = 'trimmed_data.txt',
            quote = FALSE,
            sep = '\t',
            row.names = FALSE)

