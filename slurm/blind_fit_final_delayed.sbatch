#!/bin/bash
#SBATCH --job-name=blindfit_final_delayed
#SBATCH --array=1-10
#SBATCH --ntasks=4
#SBATCH --mem=4G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=alex.s.wheaton@gmail.com

source /home/awheaton/.bash_profile
cd /home/awheaton/agn-host-project
mpirun -n 4 --oversubscribe python3 blind_fit_final.py phil_model_$SLURM_ARRAY_TASK_ID delayed_burst_final >> nohup/phil_model_$SLURM_ARRAY_TASK_ID\_final_delayed.out
