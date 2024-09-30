#!/bin/bash
#SBATCH -A m4320
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=Batch_12.log
#SBATCH --error=Batch_12.log

source ${HOME}/.bashrc

# Run the program
python3 run_batch_12.py -b 2
