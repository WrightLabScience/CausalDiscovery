#!/bin/bash

job_num=${1}
filter_var=${2}
handle_missing_var=${3}
require_trt_effect_var=${4}
search_alg_var=${5}
search_alpha_var=${6}

echo "VARIABLE VALUES:"
echo "job_num: $job_num"
echo "filter_var: $filter_var"
echo "handle_missing_var: $handle_missing_var"
echo "require_trt_effect_var: $require_trt_effect_var"
echo "search_alg_var: $search_alg_var"
echo "search_alpha_var: $search_alpha_var"
echo "" ;
echo "" ;

Rscript dataset_trimmer.R \
    --job_num $job_num \
    --filter_var $filter_var \
    --handle_missing_var $handle_missing_var \
    --require_trt_effect_var $require_trt_effect_var

output_files_prefix="Result_$job_num"

cat "variables_$job_num.txt"
cat knowledge.txt
echo "" ;

# Base tetrad command
command="java -Xmx8G -jar causal-cmd-1.12.0-jar-with-dependencies.jar  \
    --dataset trimmed_data.txt \
    --delimiter tab \
    --data-type mixed \
    --numCategories 2 \
    --knowledge knowledge.txt \
    --missing-marker NA \
    --prefix $output_files_prefix \
    --algorithm $search_alg_var"

# Add score, test, alpha, targets, etc.
# depending on algorithms used in causal search strategy
algs_need_constraint=("grasp-fci" "fci" "pc" "gfci" "rfci" "cfci" "cpc" "fci-max" "grasp" "pc-mb")
algs_need_score=("grasp-fci" "boss" "fges" "grasp" "gfci" "fges-mb")
algs_need_targets=("fges-mb" "pc-mb")

[[ " ${algs_need_constraint[@]} " =~ " $search_alg_var " ]] && command="$command --test cg-lr-test --alpha $search_alpha_var"

[[ " ${algs_need_score[@]} " =~ " $search_alg_var " ]] && command="$command --score cg-bic-score"

[[ " ${algs_need_targets[@]} " =~ " $search_alg_var " ]] && command="$command --targets time,TRT"

[ $search_alg_var == 'boss' ] && command="$command --numThreads 1"


# Print the constructed command (for debugging), then run it!
echo "Constructed command: $command"
eval $command

# remove unnecessary files
rm trimmed_data.txt
rm knowledge.txt