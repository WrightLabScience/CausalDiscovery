#!/bin/bash

java -jar ~/Desktop/EHR/causal/causal-cmd/causal-cmd-1.12.0-jar-with-dependencies.jar \
    --dataset trimmed_data.txt \
    --delimiter tab \
    --data-type mixed \
    --numCategories 2 \
    --missing-marker NA \
    --knowledge knowledge.txt \
    --algorithm grasp \
    --score cg-bic-score \
    --test cg-lr-test \
    --alpha 0.05 \
    #--score cg-bic-score \
    #--json-graph \
    #--resamplingWithReplacement \
    #--numberResampling 50