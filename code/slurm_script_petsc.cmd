#!/bin/bash

#SBATCH -o /home/hpc/pr63so/ga53lov2/petsc_console.out 
#SBATCH -D /home/hpc/pr63so/ga53lov2/proj_git/compare_solvers/code
#SBATCH -J petsc_debug
#SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=7
# node has 28 cores, and I think it has 4 sockets per node for mpp2
#SBATCH --mail-type=end 
#SBATCH --mail-user=tudor.mot@tum.de 
#SBATCH --export=NONE 
#SBATCH --time=4:00:00

./load_modules
touch /home/hpc/pr63so/ga53lov2/timing_results_petsc.txt
cd ~/proj_git/compare_solvers/code 


#this should be modified from system to system
inputfile_path=/naslx/projects/pr63so/ga53lov2/job3.sti
outputfile_path=/home/hpc/pr63so/ga53lov2/timing_results_petsc.txt

mpiexec ./time_solvers --solve-with-PETSc ${inputfile_path} ${outputfile_path} -mat_type mpisbaij -pc_type asm -ksp_type bcgs  

# -ksp_gmres_restart 200
