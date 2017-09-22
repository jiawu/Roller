#!/bin/bash     
while read P1 P2 P3 P4 P5
do
    JOB=`msub - << EOJ
#MSUB -A p20519
#MSUB -l walltime=72:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/high_sampling_2.txt
#MSUB -m bae
#MSUB -q long
#MSUB -N RandomForest
#MSUB -V        

workon seqgen
module load python/anaconda3
cd /home/jjw036/Roller/pipelines

echo "python run_tdSwing_scan.py ${P1} ${P2} ${P3} ${P4} ${P5}"
python run_tdSwing_scan.py ${P1} ${P2} ${P3} ${P4} ${P5}
EOJ
`
echo "JobID = ${JOB} for parameters ${P1} ${P2} ${P3} ${P4} ${P5} submitted on `date`"
done < job_params_high_sampling_missing.txt
exit
