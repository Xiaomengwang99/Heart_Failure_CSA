#!/bin/bash
#$ -N HF_seed
#$ -j y
#$ -t 1-64
#$ -o output_HF/$JOB_ID-$TASK_ID.log

module load R/4.3.1
Rscript HF_Survival.r
