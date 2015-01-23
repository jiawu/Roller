#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/Roller_error.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N insilico1
#MSUB -V

workon seqgen
module load python/anaconda3
cd ~/Roller
python run_pipeline_c.py
