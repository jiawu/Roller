#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=2
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/lasso_alpha_testing.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N lasso_alpha_testing
#MSUB -V

param_set=${MOAB_JOBARRAYINDEX}

## enter range 10 to get 0.01
## enter range 100 to get 0.1
## [100-1000:100]
## enter range 1000 to get 1
#param_set=$1
workon seqgen
module load python/anaconda3
cd /home/jjw036/Roller/pipelines

converted_float=${param_set}
iterating_param="alpha"
iterating_style="float ${converted_float}"

data_folder="/projects/p20519/roller_output/large_networks/Lasso/ecoli-insilico_size100_18"
output_folder="/home/jjw036/Roller/pipelines/lasso_alpha_scans/output/ecoli-insilico_size100_18_"
file_path="/home/jjw036/Roller/data/gnw_insilico/network_data/Ecoli100/Ecoli100-18_timeseries.tsv"

echo "python /home/jjw036/Roller/pipelines/run_tdSwing_scan_100_part.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}"
python /home/jjw036/Roller/pipelines/run_tdSwing_scan_100_part.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}
#python run_tdSwing_scan_vgranger.py ${data_folder} ${output_folder} ${file_path} ${iterating_param2} ${iterating_style2}

#data_folder=${data_folder/Dionesus/Dionesus}
#output_folder=${output_folder/Dionesus/Dionesus}
#python run_tdDionesus_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}

#data_folder=${data_folder/Dionesus/Dionesus}
#output_folder=${output_folder/Dionesus/Dionesus}
#python run_tdDionesus_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}

#python run_tdDionesus_scan.py /projects/p20519/roller_output/optimizing_window_size/Swing/janes /projects/p20519/roller_output/stability_analysis/Dionesus/janes_${iterating_param}_ /projects/p20519/Roller/data/invitro/janes_timeseries.tsv n_trees log
#python run_pipeline_RF_window_scan_janes.py ${nwindows}
