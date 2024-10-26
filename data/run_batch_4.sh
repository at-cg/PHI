#!/bin/bash
#SBATCH -A m4320
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --output=Batch_4.log
#SBATCH --error=Batch_4.log

source ${HOME}/.bashrc

# Run the program
python3 run_batch_4.py -b 2
