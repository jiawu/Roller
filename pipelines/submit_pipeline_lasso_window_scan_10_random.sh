#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/Roller_error.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N lasso_window_scan_10
#MSUB -V

nwindows=${MOAB_JOBARRAYINDEX}

workon seqgen
module load python/anaconda3
cd ~/Roller
python run_pipeline_lasso_window_scan_10_random.py ${nwindows}
