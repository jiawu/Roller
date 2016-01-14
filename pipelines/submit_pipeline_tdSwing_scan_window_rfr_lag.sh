#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/td_scan.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N lag_scans for dream4
#MSUB -V

param_set=${MOAB_JOBARRAYINDEX}
#param_set=$1
workon seqgen
module load python/anaconda
cd /home/jjw036/Roller/pipelines


iterating_param="max_lag"
iterating_style="num 1 2 3 4 5"

iterating_param2="min_lag"
iterating_style2="num 0 1 2 3"

if [[ ${param_set} == 1 ]] 
then
    #echo "param_set is ${param_set}"
    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/janes"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/janes_${iterating_param}_"
    file_path="/projects/p20519/Roller/data/invitro/janes_timeseries.tsv"
elif [[ ${param_set} == 2 ]]
then
    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/whitfield_muk"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/whitfield_muk_${iterating_param}_"
    file_path="/projects/p20519/Roller/data/invitro/whitfield_muk_timeseries.tsv"
elif [[ ${param_set} == 3 ]]
then
    data_folder="/projects/p20519/roller_output/optimizing_window_size/RandomForest/whitfield_shojaie"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/whitfield_shojaie_${iterating_param}_"
    file_path="/projects/p20519/Roller/data/invitro/whitfield_shojaie_timeseries.tsv"
elif [[ "${param_set}" -ge 4 && "${param_set}" -lt 9 ]]
then
    #process the in silico datasets a little bit more smoothly
    #param_set 4 = dataset 1, param_set 5 = dataset 2, etc etc
    insilico_dataset_index=$((${param_set}-3))

    data_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/insilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/stability_analysis/RandomForest/insilico_size10_${insilico_dataset_index}_ntrees_"
    file_path="/home/jjw036/Roller/data/dream4/insilico_size10_${insilico_dataset_index}_timeseries.tsv"
    
elif [[ "${param_set}" -ge 9 && "${param_set}" -lt 14 ]]
then
    insilico_dataset_index=$((${param_set}-8))

    data_folder="/projects/p20519/roller_output/sampling/RandomForest/uniform_samplinginsilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/sampling/RandomForest/uniform_samplinginsilico_size10_${insilico_dataset_index}_ntrees_"
    file_path="/home/jjw036/Roller/data/dream4/uniform_sampling/uniform_samplinginsilico_size10_${insilico_dataset_index}_timeseries.tsv"

elif [[ ${param_set} -ge 14 && ${param_set} -lt 19 ]]
then
    insilico_dataset_index=$((${param_set}-13))

    data_folder="/projects/p20519/roller_output/sampling/RandomForest/nonuniform_samplinginsilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/sampling/RandomForest/nonuniform_samplinginsilico_size10_${insilico_dataset_index}_ntrees_"
    file_path="/home/jjw036/Roller/data/dream4/nonuniform_sampling/nonuniform_samplinginsilico_size10_${insilico_dataset_index}_timeseries.tsv"

elif [[ ${param_set} -ge 19 && ${param_set} -lt 39 ]]
then
    insilico_dataset_index=$((${param_set}-18))

    data_folder="/projects/p20519/roller_output/gnw/RandomForest/yeast-insilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/gnw/RandomForest/yeast-insilico_size10_${insilico_dataset_index}_"
    file_path="/home/jjw036/Roller/data/gnw_insilico/network_data/Yeast/Yeast-${insilico_dataset_index}_timeseries.tsv"

elif [[ ${param_set} -ge 39 && ${param_set} -lt 59 ]]
then
    insilico_dataset_index=$((${param_set}-38))

    data_folder="/projects/p20519/roller_output/gnw/RandomForest/ecoli-insilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/gnw/RandomForest/ecoli-insilico_size10_${insilico_dataset_index}_"
    file_path="/home/jjw036/Roller/data/gnw_insilico/network_data/Ecoli/Ecoli-${insilico_dataset_index}_timeseries.tsv"

elif [[ ${param_set} -ge 59 && ${param_set} -lt 60 ]]
then

    data_folder="/projects/p20519/roller_output/dream8/RandomForest/insilico"
    output_folder="/projects/p20519/roller_output/dream8/RandomForest/insilico"
    file_path="/home/jjw036/Roller/data/dream8/insilico/insilico_timeseries.tsv"

elif [[ ${param_set} -ge 60 && ${param_set} -lt 92 ]]
then
    insilico_dataset_index=$((${param_set}-60))
    ind=$((${insilico_dataset_index}%8))

    if [[ ${ind} -eq 0 ]]
    then
        perturb="EGF"
    elif [[ ${ind} -eq 1 ]]
    then
        perturb="FGF"    
    elif [[ ${ind} -eq 2 ]]
    then
        perturb="HGF"
    elif [[ ${ind} -eq 3 ]]
    then
        perturb="IGF1"
    elif [[ ${ind} -eq 4 ]]
    then
        perturb="Insulin"
    elif [[ ${ind} -eq 5 ]]
    then
        perturb="NRG1"
    elif [[ ${ind} -eq 6 ]]
    then
        perturb="PBS"
    elif [[ ${ind} -eq 7 ]]
    then
        perturb="Serum"
    fi

    if [[ "${insilico_dataset_index}" -ge 0 && "${insilico_dataset_index}" -lt 8 ]]
    then
        cell_type="UACC812"
    elif [[ ${insilico_dataset_index} -ge 8 && "${insilico_dataset_index}" -lt 16 ]]
    then
        cell_type="MCF7"
    elif [[ ${insilico_dataset_index} -ge 16 && "${insilico_dataset_index}" -lt 24 ]]
    then
        cell_type="BT549"
    elif [[ ${insilico_dataset_index} -ge 24 && "${insilico_dataset_index}" -lt 32 ]]
    then
        cell_type="BT20"
    fi
    data_folder="/projects/p20519/roller_output/dream8/RandomForest/${cell_type}_${perturb}_"
    output_folder="/projects/p20519/roller_output/dream8/RandomForest/${cell_type}_${perturb}_"
    file_path="/home/jjw036/Roller/data/dream8/invitro/${cell_type}_${perturb}_timeseries.tsv"

elif [[ ${param_set} -ge 92 && ${param_set} -lt 112 ]]
then
    insilico_dataset_index=$((${param_set}-91))

    data_folder="/projects/p20519/roller_output/community/yeast-insilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/community/yeast-insilico_size10_${insilico_dataset_index}_"
    file_path="/home/jjw036/Roller/data/gnw_insilico/network_data/Yeast/Yeast-${insilico_dataset_index}_timeseries.tsv"

elif [[ ${param_set} -ge 112 && ${param_set} -lt 132 ]]
then
    insilico_dataset_index=$((${param_set}-111))

    data_folder="/projects/p20519/roller_output/community/ecoli-insilico_size10_${insilico_dataset_index}"
    output_folder="/projects/p20519/roller_output/community/ecoli-insilico_size10_${insilico_dataset_index}_"
    file_path="/home/jjw036/Roller/data/gnw_insilico/network_data/Ecoli/Ecoli-${insilico_dataset_index}_timeseries.tsv"
  

  
else
    echo "out of bounds"
fi

echo "python run_tdSwing_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}"
python run_tdSwing_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}
python run_tdSwing_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param2} ${iterating_style2}

#data_folder=${data_folder/RandomForest/Lasso}
#output_folder=${output_folder/RandomForest/Lasso}
#python run_tdRandomForest_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}

#data_folder=${data_folder/RandomForest/Dionesus}
#output_folder=${output_folder/RandomForest/Dionesus}
#python run_tdRandomForest_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}

#python run_tdRandomForest_scan.py /projects/p20519/roller_output/optimizing_window_size/Swing/janes /projects/p20519/roller_output/stability_analysis/Lasso/janes_${iterating_param}_ /projects/p20519/Roller/data/invitro/janes_timeseries.tsv n_trees log
#python run_pipeline_RF_window_scan_janes.py ${nwindows}
