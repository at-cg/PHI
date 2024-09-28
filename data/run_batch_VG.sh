#!/bin/bash
#SBATCH -A m4320
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --output=Batch_VG.log
#SBATCH --error=Batch_VG.log

source ${HOME}/.bashrc

# Run the program
python3 run_VG.py -b 4
