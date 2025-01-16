#!/bin/bash

job_num=${1}
filter_var=${2}
select_var=${3}
handle_missing_var=${4}
include_outcome_var=${5}
search_alg_var=${6}
search_alpha_var=${7}

echo "VARIABLE VALUES:"
echo "job_num: $job_num"
echo "filter_var: $filter_var"
echo "select_var: $select_var"
echo "handle_missing_var: $handle_missing_var"
echo "include_outcome_var: $include_outcome_var"
echo "search_alg_var: $search_alg_var"
echo "search_alpha_var: $search_alpha_var"
echo "" ;
echo "" ;

Rscript dataset_trimmer.R \
    --filter_var $filter_var \
    --select_var $select_var \
    --handle_missing_var $handle_missing_var \
    --include_outcome_var $include_outcome_var

output_files_prefix="Result_$job_num"

cat knowledge.txt
echo "" ;

java -Xmx8G -jar causal-cmd-1.12.0-jar-with-dependencies.jar \
    --dataset trimmed_data.txt \
    --delimiter tab \
    --data-type mixed \
    --numCategories 2 \
    --knowledge knowledge.txt \
    --algorithm $search_alg_var \
    --alpha $search_alpha_var \
    --test cg-lr-test \
    --prefix $output_files_prefix # --json-graph

rm trimmed_data.txt
rm knowledge.txt