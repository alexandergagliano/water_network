#!/bin/tcsh
##### These lines are for Moab
#  predict the duration of the job
#SBATCH -t 24:00:00  # Walltime
#SBATCH -N 1       # Number of nodes
#SBATCH -o slurm_%j.out # name of the stdout
#SBATCH --qos=standby
#SBATCH -J PopIb       # job name
#SBATCH -d singleton    #ensures one job starts after another
#SBATCH -A w18_lgsims
#
#  specify the pathname for output
#
#  combine stdout and stdin
##MSUB -j oe
#
#  forward current environment variables to the job
#MSUB -V
date
#module purge
#module load gcc/5.3.0 
#module load openmpi/1.10.5

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/lustre/scratch3/turquoise/agagliano/WATER/lib
python production_freefall_A_plotAbundances.py
date
