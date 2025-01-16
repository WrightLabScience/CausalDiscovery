#!/bin/bash

java -jar ~/Desktop/EHR/causal/causal-cmd/causal-cmd-1.12.0-jar-with-dependencies.jar \
    --dataset datasets/data.txt \
    --delimiter tab \
    --data-type mixed \
    --numCategories 2 \
    --knowledge knowledge.txt \
    --algorithm pc \
    --alpha 0.05 \
    --test cg-lr-test \
    --json-graph \
    --out runs
    #--score cg-bic-score \
    #--json-graph \
    #--resamplingWithReplacement \
    #--numberResampling 50