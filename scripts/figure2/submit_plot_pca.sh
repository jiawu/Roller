#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=4
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/pca.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N PCA
#MSUB -V

workon seqgen
module load python/anaconda3
cd /home/jjw036/Roller/scripts/figure2/

python plot_fig_2b.py
