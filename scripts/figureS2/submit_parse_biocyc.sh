#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/insilico_rank_plots.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N parse_biocyc_all
#MSUB -V

param_set=${MOAB_JOBARRAYINDEX}
#param_set=$1
workon seqgen
module load python/anaconda3
source activate my_root
cd /home/jjw036/Roller/scripts/figure4

python parse_biocyc.py community 26
