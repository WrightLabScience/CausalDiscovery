# pull tetrad search results from the grid after running grid search
# do some analysis here

setwd('~/Desktop/tetrad_cmd_files/MRSAbact/')
lf <- list.files(path = 'Results/', pattern='txt')

params <- read.table('job_files/search.map',
                     sep=' ', 
                     row.names=1,
                     col.names = c('key', 'filter_var', 'select_var', 'handle_missing_var', 'include_outcome_var', 'search_alg_var', 'search_alpha_var'))

# question number 0: was the temporal order followed? does time cause anything? does TRT cause anything except time?
# question number 1: how many and which runs revealed causes of TRT and TRT causing time?
# distribution of edges types?

for (i in seq_along(lf)) {
   print(params[i,])
   cat('\n')
   
   out <- readLines(con = paste0('Results/', lf[i]))
   edges <- out[(grep('Graph Edges:', out) + 1):(grep('Graph Attributes:', out) - 2)]
   edges <- gsub('^[0-9]+\\. ', '', edges)
   edge_types <- gsub('[A-z0-9]+ (.+) [A-z0-9]+', '\\1', edges)
   edge_list <- strsplit(edges, split=' [-o<>]{3} ')
   
   cat('Total number of edges:', length(edge_list), '\n')
   # distribution of edge types
   cat('Distribution of ')
   print(table(edge_types))
   cat('\n')
   
   
   # check for temporal order violations?
   if (any(sapply(edge_list, function(x) x[1] == 'time')))
      cat('Warning: time causes something\n')
   if (any(sapply(edge_list, function(x) x[1] == 'TRT' & x[2] != 'time')))
      cat('Warning: TRT causes something other than time\n')
   
   
   # does TRT cause time?
   if (!any(sapply(edge_list, function(x) x[1] == 'TRT' & x[2] == 'time'))) {
      cat('Warning: TRT does not cause time\n') 
   } else {
      cat('TRT causes time!!\n') 
   }
   
   
   # what causes TRT and time or both?
   cause_trt <- cause_surv <- character()
   for (pair in edge_list) {
      if (pair[1] != 'time' & pair[2] == 'TRT')
         cause_trt <- c(cause_trt, pair[1])
      if (pair[1] != 'TRT' & pair[2] == 'time')
         cause_surv <- c(cause_surv, pair[2])
   }
   if (length(cause_trt) > 0) {
      cat('CAUSE TRT:', cause_trt, '\n')
   } else cat('Warning: no causes of TRT!\n')
   if (length(cause_surv) > 0) {
      cat('CAUSE SURV:', cause_surv, '\n')
   } else cat('Warning: no causes of time!\n')
   cause_conf <- intersect(cause_trt, cause_surv)
   if (length(cause_conf) > 0) {
      cat('CAUSE BOTH:', cause_conf, '\n')
   } else cat('Warning: no common causes!\n')
   
   
   cat('\n__________________________________________\n\n')
}























