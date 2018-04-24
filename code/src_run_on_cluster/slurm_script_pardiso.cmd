#!/bin/bash
 

#SBATCH -o /home/hpc/pr63so/ga53lov2/pard_debug_inf_err.%j.%N.out 
#SBATCH -D /home/hpc/pr63so/ga53lov2/code
#SBATCH -J pard_8thr_quicktest
#SBATCH --get-user-env 
#SBATCH --clusters=mpp2
#SBATCH --nodes=1-1 
#SBATCH --cpus-per-task=28 
#SBATCH --mail-type=end 
#SBATCH --mail-user=tudor.mot@tum.de 
#SBATCH --export=NONE 
#SBATCH --time=1:00:00 

./load_modules
touch outputfile_path=/home/hpc/pr63so/ga53lov2/timing_results_pardiso.txt
cd ~/code
export OMP_NUM_THREADS=28
# 28 is the maximum reasonable value for CooLMUC-2 

#this should be modified from system to system
inputfile_path=/naslx/projects/pr63so/ga53lov2/job3.sti
outputfile_path=/home/hpc/pr63so/ga53lov2/timing_results_pardiso.txt

mpiexec -n 1 ./release --solve-with-pardiso ${inputfile_path} ${outputfile_path} -ksp_type gmres -ksp_gmres_restart 200