#!/bin/bash
#SBATCH --job-name=blindfit_wide_model_06
#SBATCH --ntasks=4
#SBATCH --mem=8G
#SBATCH --time=36:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=alex.s.wheaton@gmail.com

source /home/awheaton/.bash_profile
cd /home/awheaton/agn-host-project
mpirun -n 4 --oversubscribe nohup python3 blind_fit_wide.py phil_model_06 > nohup/phil_model_06_wide.out
