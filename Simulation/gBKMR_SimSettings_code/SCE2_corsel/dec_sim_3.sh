#!/bin/bash
#$ -cwd -l mem=3g,time=18:10:00 -S /bin/bash -N JOBs3 -j y -t 1-500

currind=$SGE_TASK_ID

R=/nfs/apps/R/3.6.0/bin/R
R_LIBS_USER=/ifs/home/msph/LeeLab/zc2326/R_LIB:/ifs/scratch/msph/software/R/library360:/ifs/scratch/msph/software/R/library:$R_LIBS_USER

${R} --vanilla --args $currind < dec_sim_3.R 
