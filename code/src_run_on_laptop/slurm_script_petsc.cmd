#!/bin/bash

#SBATCH -o /home/hpc/pr63so/ga53lov2/petsc_small_debug_2proc.%j.%N.out 
#SBATCH -D /home/hpc/pr63so/ga53lov2/code
#SBATCH -J petsc_debug
#SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --ntasks=28
# multiples of 28 for mpp2
#SBATCH --mail-type=end 
#SBATCH --mail-user=tudor.mot@tum.de 
#SBATCH --export=NONE 
#SBATCH --time=8:00:00 

./load_modules
touch /home/hpc/pr63so/ga53lov2/timing_results_petsc.txt
cd ~/code 

#this should be modified from system to system
inputfile_path=/naslx/projects/pr63so/ga53lov2/beamp.txt
outputfile_path=/home/hpc/pr63so/ga53lov2/timing_results_petsc.txt

mpiexec ./release --solve-with-PETSc ${inputfile_path} ${outputfile_path} -ksp_type bcgs -ksp_monitor -info -info_exclude pc,ksp 

# -ksp_gmres_restart 200
