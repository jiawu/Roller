#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/insilico_rank_plots.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N lasso_alpha
#MSUB -V

param_set=${MOAB_JOBARRAYINDEX}
#param_set=$1
workon seqgen
module load python/anaconda3
source activate my_root
cd /home/jjw036/Roller/scripts/figureS2

python lasso_alpha_selection_10b.py Ecoli/Ecoli-${param_set}
