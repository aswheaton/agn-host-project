#!/bin/bash
#SBATCH --job-name=blindfit_wide
#SBATCH --array=01,02,03,04,05,06,07,08,09,10
#SBATCH --ntasks=4
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=alex.s.wheaton@gmail.com

source /home/awheaton/.bash_profile
cd /home/awheaton/agn-host-project
mpirun -n 4 --oversubscribe nohup python3 blind_fit_wide.py phil_model_$SLURM_ARRAY_TASK_ID > nohup/phil_model_$SLURM_ARRAY_TASK_ID\_.wide.out
