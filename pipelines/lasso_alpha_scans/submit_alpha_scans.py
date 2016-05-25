#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=72:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/alpha.txt
#MSUB -m bae
#MSUB -q long
#MSUB -N alpha_select
#MSUB -V

workon seqgen
module load python/anaconda3
cd /home/jjw036/Roller/pipelines/lasso_alpha_scans/

python alpha_scans.py

