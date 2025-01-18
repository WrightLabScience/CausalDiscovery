# pull tetrad search results from the grid after running grid search
# do some analysis here

setwd('~/Desktop/CausalDiscovery/')
lf <- list.files(path = 'Results/', pattern='txt')

params <- read.table('job_files/search.map',
                     sep=' ', 
                     row.names=1,
                     col.names = c('key', 'filter_var', 'select_var', 'handle_missing_var', 'include_outcome_var', 'search_alg_var', 'search_alpha_var'))

# question number 0: was the temporal order followed? does time cause anything? does TRT cause anything except time?
# question number 1: how many and which runs revealed causes of TRT and TRT causing time?
# distribution of edges types?

calcVennCounts <- function(df, edge_type) {
   cat('Counts for edge_type: ', edge_type, '\n')
   
   cause_trt <- df$cause[df$cause != 'time' & df$effect == 'TRT' & df$type == edge_type]
   cause_time <- df$cause[df$cause != 'TRT' & df$effect == 'time' & df$type == edge_type]
   cause_both <- intersect(cause_trt, cause_time)
   
   cat('Num causes of TRT:\t', length(cause_trt), '\n')
   cat('Num causes of time:\t', length(cause_time), '\n')
   cat('Num confounders:\t', length(cause_both), ifelse(length(cause_both) > 0L, paste(cause_both, collapse=' '), ''), '\n')
   cat('\n')
}

for (i in seq_along(lf)) {
   print(params[i,])
   cat('\n')
   
   out <- readLines(con = paste0('Results/', lf[i]))
   start_edges <- grep('Graph Edges:', out) + 1
   end_edges <- grep('Graph Attributes:', out) - 2
   if (length(end_edges) == 0L)
      end_edges <- length(out)
   edges <- out[start_edges:end_edges]
   edges <- gsub('^[0-9]+\\. ', '', edges)
   edge_types <- gsub('[A-z0-9]+ (.+) [A-z0-9]+', '\\1', edges)
   edge_list <- strsplit(edges, split=' [-o<>]{3} ')
   df <- data.frame(
      cause = sapply(edge_list, '[', 1),
      effect = sapply(edge_list, '[', 2),
      type = edge_types
   )
   all_variables <- unique(c(df$cause, df$effect))
   
   cat('Total number of edges:', nrow(df), '\n')
   # distribution of edge types
   cat('Distribution of edge types:')
   print(table(df$type))
   cat('\n')
   
   
   # check for temporal order violations?
   if (any(df$cause == 'time' | (df$cause == 'TRT' & df$effect != 'time'))) {
      cat('Error: temporal violation.\n')
   }
   
   # does TRT cause time?
   if (any(df$cause == 'TRT' & df$effect == 'time')) {
      cat('Success: TRT causes time.\n')
   } else {
      cat('Warning: TRT does not cause time.\n') 
   }
   
   for (edge_type in unique(df$type)) {
      calcVennCounts(df, edge_type)
   }
   
   cat('\n__________________________________________\n\n')
}























