#!/bin/tcsh
##### These lines are for Moab
#  predict the duration of the job
#SBATCH -t 00:10:00  # Walltime
#SBATCH -N 1       # Number of nodes
#SBATCH -o slurm_%j.out # name of the stdout
#SBATCH -J Pop6       # job name
##SBATCH -d singleton    #ensures one job starts after another
#SBATCH -A w18_lgsims
#
#  specify the pathname for output
#
#  combine stdout and stdin
##MSUB -j oe
#
#  forward current environment variables to the job
##MSUB -V
date
#module purge
#module load gcc/5.3.0 
#module load openmpi/2.1.2
#module load python/3.5-anaconda-4.1.1

#setenv LD_LIBRARY_PATH /lustre/scratch3/turquoise/agagliano/WATER/lib:${LD_LIBRARY_PATH}
#setenv PYTHONPATH /lustre/scratch3/turquoise/agagliano/WATER/Fresh_1027/src/python/:${PYTHONPATH}
#setenv PATH /users/agagliano/.conda/envs/my_root/bin/python:/lustre/scratch3/turquoise/agagliano/WATER/yt-conda/bin:${PATH}

python oneZone_UMIST2012.py
date
