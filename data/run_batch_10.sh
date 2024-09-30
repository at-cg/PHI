#!/bin/bash
#SBATCH -A m4320
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=Batch_10.log
#SBATCH --error=Batch_10.log

source ${HOME}/.bashrc

# Run the program
python3 run_batch_10.py -b 2
