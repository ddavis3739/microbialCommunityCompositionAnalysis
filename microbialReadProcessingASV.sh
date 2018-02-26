#!/bin/bash

## STEPS, following (steps from microbialReadProcessingOTU.sh)
# 7) create deblur tables

# run deblur workflow, creates table with representative sequences instead of OTU label
# --min-reads default is 10 but is 2 to match with the OTU picking method done in this workflow
deblur workflow --seqs-fp combined/seqs.fna --output-dir deblur -t 300 --min-reads 2

# move output to accesible location
cp deblur/reference-hit.biom ../reference-hit_all.biom