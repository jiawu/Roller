#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=140:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/td_scan.txt
#MSUB -m bae
#MSUB -q long
#MSUB -N td_scan_ntrees_log
#MSUB -V

param_set=${MOAB_JOBARRAYINDEX}
#param_set=$1
workon seqgen
module load python/anaconda
cd /home/jjw036/Roller/pipelines


iterating_param="sort_by"
iterating_style="string mean rank adj"
#iterating_style="string min_min min_max min_mean min_median max_min max_max max_mean max_median mean_min mean_max mean_mean mean_median median_min median_max median_mean mean_median median_min median_max median_mean median_median"
if [ $param_set -eq 1 ] 
then
    #echo "param_set is ${param_set}"
    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/janes"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/janes_${iterating_param}_"
    file_path="/projects/p20519/Roller/data/invitro/janes_timeseries.tsv"
elif [ ${param_set} == 2 ]
then
    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/whitfield_muk"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/whitfield_muk_${iterating_param}_"
    file_path="/projects/p20519/Roller/data/invitro/whitfield_muk_timeseries.tsv"
elif [ ${param_set} == 3 ]
then
    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/whitfield_shojaie"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/whitfield_shojaie_${iterating_param}_"
    file_path="/projects/p20519/Roller/data/invitro/whitfield_shojaie_timeseries.tsv"
elif [ ${param_set} -ge 4 ]
then
    #process the in silico datasets a little bit more smoothly
    #param_set 4 = dataset 1, param_set 5 = dataset 2, etc etc
    insilico_dataset_index=$((${param_set}-3))

    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/insilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/insilico_size10_${insilico_dataset_index}_ntrees_"
    file_path="/projects/p20519/Roller/data/dream4/insilico_size10_${insilico_dataset_index}_timeseries.tsv"
    
fi

echo "python run_tdSwing_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}"
python run_tdSwing_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}
#python run_tdSwing_scan.py /projects/p20519/roller_output/optimizing_window_size/RandomForest/janes /projects/p20519/roller_output/stability_analysis/RandomForest/janes_${iterating_param}_ /projects/p20519/Roller/data/invitro/janes_timeseries.tsv n_trees log
#python run_pipeline_RF_window_scan_janes.py ${nwindows}

