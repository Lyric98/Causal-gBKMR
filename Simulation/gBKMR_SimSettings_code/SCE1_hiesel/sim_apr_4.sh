#!/bin/bash
#$ -cwd -l mem=12g,time=18:10:00 -S /bin/bash -N JOBs4 -j y  

currind=$SGE_TASK_ID

R=/nfs/apps/R/3.6.0/bin/R
export R_LIBS_USER=/nfs/apps/R/3.6.0/lib64/R/library:/ifs/scratch/msph/software/R/library351:/ifs/scratch/msph/software/R/library351_Exp:/ifs/home/msph/biostat/zc2326/R_LIB:/ifs/scratch/msph/software/R/library360:/ifs/scratch/msph/software/R/library:$R_LIBS_USER

${R} --vanilla --args $currind < stantest.R 
