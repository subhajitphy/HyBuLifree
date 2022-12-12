#!/bin/sh
#PBS -N job_ew
#PBS -o op_ew.txt
#PBS -e er_ew.txt
#PBS -l select=1:ncpus=32:mem=72gb
#PBS -q small
echo "Starting at "`date`
source /scratch/subhajit.dandapat/Packages/miniconda3/bin/activate
. /scratch/subhajit.dandapat/Packages/lalsuite/_inst/etc/lalsuite-user-env.sh
conda activate gwhyp
cd $PBS_O_WORKDIR
export PATH=$PATH:/home/apps/GnuParallel/bin
echo $PBS_JOBID
python run_cons_e.py
echo "finishing at "`date`
