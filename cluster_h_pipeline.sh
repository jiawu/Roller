#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=72:00:00
#MSUB -l nodes=1:ppn=2
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/Roller_error_h.txt
#MSUB -m bae
#MSUB -q long
#MSUB -N cluster_h
#MSUB -V

nwindows=${MOAB_JOBARRAYINDEX}

workon seqgen
module load python/anaconda
cd ~/Roller
python run_pipeline_h.py ${nwindows}
