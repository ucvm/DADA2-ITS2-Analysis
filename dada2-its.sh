#! /usr/bin/env bash

#BSUB -J dada2-its2
#BSUB -n 4 
#BSUB -R "span[hosts=1]"
#BSUB -We 20:00
#BSUB -o test_%J-%I.out
#BSUB -e test_%J-%I.err

source activate qiime2-2018.11
Rscript dada2-its.r
source deactivate
