#!/bin/bash

#SBATCH --account=biostats
#SBATCH -c 1 # Number of cores
#SBATCH -N 1
#SBATCH -t 3-23:00
#SBATCH -o /burg/biostats/users/yl5465/Causal-gBKMR/Simulation/cluster/s6_result/test_error.out #specify where to save errors returned by the program
#SBATCH -e /burg/biostats/users/yl5465/Causal-gBKMR/Simulation/cluster/s6_result/test_log.err #specify where to save the output log
#SBATCH --array=1-100 #number of jobs to run, it is currently set to 1 job(1 year), change it to array=1-13 for 13 years(jobs)
#SBATCH --mem=4g #memory requested
#SBATCH -J gBKMR  #job name, this case:aw-area weighted aggregation
#SBATCH --mail-type=ALL #notifications for job done
#SBATCH --mail-user=yl5465@cumc.columbia.edu # send to address

module load R 


R CMD BATCH --quiet --no-restore --no-save s6_corsel.R /burg/biostats/users/yl5465/Causal-gBKMR/Simulation/cluster/s6_result/s6_${SLURM_ARRAY_TASK_ID}.Rout
