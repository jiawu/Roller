#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=2
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/omranian_dionesus.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N dionesus_largenetworks
#MSUB -V

param_set=${MOAB_JOBARRAYINDEX}
#param_set=$1
workon seqgen
module load python/anaconda3
cd /home/jjw036/Roller/pipelines


iterating_param="td_window"
iterating_style="num ${param_set}"

data_folder="/projects/p20519/roller_output/large_networks/Dionesus/omranian"
output_folder="/projects/p20519/roller_output/large_networks/Dionesus/omranian_"
file_path="/home/jjw036/Roller/data/invitro/omranian_parsed_timeseries.tsv"

echo "python run_tdSwing_scan_custom.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}"
python run_tdSwing_scan_custom.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}
#python run_tdSwing_scan_vgranger.py ${data_folder} ${output_folder} ${file_path} ${iterating_param2} ${iterating_style2}

#data_folder=${data_folder/Dionesus/Dionesus}
#output_folder=${output_folder/Dionesus/Dionesus}
#python run_tdDionesus_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}

#data_folder=${data_folder/Dionesus/Dionesus}
#output_folder=${output_folder/Dionesus/Dionesus}
#python run_tdDionesus_scan.py ${data_folder} ${output_folder} ${file_path} ${iterating_param} ${iterating_style}

#python run_tdDionesus_scan.py /projects/p20519/roller_output/optimizing_window_size/Swing/janes /projects/p20519/roller_output/stability_analysis/Dionesus/janes_${iterating_param}_ /projects/p20519/Roller/data/invitro/janes_timeseries.tsv n_trees log
#python run_pipeline_RF_window_scan_janes.py ${nwindows}
