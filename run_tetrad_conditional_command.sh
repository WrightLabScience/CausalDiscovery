#!/bin/bash


# Base command
command="java -jar ~/Desktop/EHR/causal/causal-cmd/causal-cmd-1.12.0-jar-with-dependencies.jar  \
    --dataset datasets/data.txt \
    --delimiter tab \
    --data-type mixed \
    --numCategories 2 \
    --knowledge knowledge.txt \
    --json-graph \
    --out runs \
    --prefix $output_files_prefix \
    --algorithm $search_alg_var"

# Add score, test, and alpha, depending on algorithms used in causal search strategy
algs_need_constraint=("grasp-fci" "fci" "pc" "gfci" "rfci" "cpc" "fci-max")
algs_need_score=("grasp-fci" "boss" "fges" "grasp" "gfci")

if [[ " ${algs_need_constraint[@]} " =~ " $search_alg_var " ]]; then
    echo "alg needs a test "
    command="$command --test cg-lr-test --alpha $search_alpha_var"
fi

if [[ " ${algs_need_score[@]} " =~ " $search_alg_var " ]]; then
    echo "alg needs a score "
    command="$command --score cg-bic-score"
fi

# Print the constructed command (for debugging)
echo "Constructed command: $command"

# Execute the causal-cmd command
eval $command

rm trimmed_data.txt