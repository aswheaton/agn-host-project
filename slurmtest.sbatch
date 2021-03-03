#!/bin/bash
#SBATCH --job-name=slurmtest
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCH --time=00:01:00
#SBATCH --mail-type=END
#SBATCH --mail-user=alex.s.wheaton@gmail.com

source /home/awheaton/.bash_profile
cd /home/awheaton/agn-host-project
nohup python3 slurmtest.py > nohup/slurmtest.out
