#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=99:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/Roller_error.txt
#MSUB -m bae
#MSUB -q long
#MSUB -N RF_window_scan_100
#MSUB -V

nwindows=${MOAB_JOBARRAYINDEX}

workon seqgen
module load python/anaconda
cd ~/Roller
python run_pipeline_RF_window_scan_100.py ${nwindows}
