#!/bin/bash
#SBATCH -A m4320
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=Batch_11.log
#SBATCH --error=Batch_11.log

source ${HOME}/.bashrc

# Run the program
python3 run_batch_11.py -b 2
