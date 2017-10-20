#!/bin/bash     
while read P1 P2 P3 P4 P5
do
    JOB=`msub - << EOJ
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/large.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N fine_sampling
#MSUB -V        

workon seqgen
module load python/anaconda3
cd /home/jjw036/Roller/pipelines

echo "python run_tdSwing_scan.py ${P1} ${P2} ${P3} ${P4} ${P5}"
python run_tdSwing_scan.py ${P1} ${P2} ${P3} ${P4} ${P5}
EOJ
`
echo "JobID = ${JOB} for parameters ${P1} ${P2} ${P3} ${P4} ${P5} submitted on `date`"
done < job_params_min_max_tdwindow.txt
exit
