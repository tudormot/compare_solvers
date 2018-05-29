#!/bin/bash
 

#SBATCH -o /home/hpc/pr63so/ga53lov2/pard_console.out 
#SBATCH -D /home/hpc/pr63so/ga53lov2/proj_git/compare_solvers/code
#SBATCH -J pard_8thr_quicktest
#SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --nodes=1-1 
#SBATCH --cpus-per-task=4 
#SBATCH --mail-type=end 
#SBATCH --mail-user=tudor.mot@tum.de 
#SBATCH --export=NONE 
#SBATCH --time=1:00:00 

./load_modules
outputfile_path=/home/hpc/pr63so/ga53lov2/timing_results_pardiso.txt
touch ${outputfile_path}

cd ~/proj_git/compare_solvers/code
export OMP_NUM_THREADS=4
# 28 is the maximum reasonable value for CooLMUC-2 

#this should be modified from system to system
inputfile_path=/naslx/projects/pr63so/ga53lov2/job3.sti
outputfile_path=/home/hpc/pr63so/ga53lov2/timing_results_pardiso.txt

mpiexec -n 1 ./time_solvers --solve-with-pardiso ${inputfile_path} ${outputfile_path} 
