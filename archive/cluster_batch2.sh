#!/bin/bash
#MSUB -A p20519
#MSUB -l walltime=24:00:00
#MSUB -l nodes=1:ppn=1
#MSUB -M jiawu@u.northwestern.edu
#MSUB -j oe
#MSUB -o /projects/p20519/jia_output/error.txt
#MSUB -m bae
#MSUB -q normal
#MSUB -N insilico1
#MSUB -V

export R_LIBS="/home/jjw036/R/library"
export PATH="$PATH:/home/jjw036/.local/bin"
#:/home/jjw036/meme/bin:
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/jjw036/R/library"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/jjw036/.local/lib"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/jjw036/.local/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/bin"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/hpc/software"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/home/jjw036/.local/lib64/R/lib"
export WORKON_HOME="/home/jjw036/virtualenv"
export VIRTUALENVWRAPPER_PYTHON="/software/enthought/python/epd_free-7.3.2-rh6-x86_64/bin/python2.7"
export VIRTUALENVWRAPPER_VIRTUALENV="/home/jjw036/.local/bin/virtualenv"
. /home/jjw036/.local/bin/virtualenvwrapper.sh

workon seqgen
module load python/anaconda3
nbatch=${MOAB_JOBARRAYINDEX}

python run_pipeline_a.py -i dream_2_params.txt -w ${MOAB_JOBARRAYINDEX} -s insilico2_${MOAB_JOBARRAYINDEX}.result
