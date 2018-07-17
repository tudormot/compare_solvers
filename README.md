# compare_solvers
1. Short description:
  The aim of this repository is to provide tools to easily compare the performance of linear system solver libraries. For the time being only two such libraries are considered: intel MKL (with its direct solver Intel Pardiso) and PETSc (any solver).
  The repository contains the following:
  1.1 The "compare_solvers" code, IE the main program, which can be found in directory /src
  1.2 Python script to parse the "compare_solvers" output, and generate some nice graphs for solver comparison. The script is available in /util/python_scripts. The script has it's own README, again available in /util/python_scripts
  1.3 A mini-program to compare solutions once they were stored in files, available in /util/sol_compare. This was mostly useful during debugging/development. The main program, "compare_solvers", already has a compare solutions built-in functionality, and it is much easier to use that rather than saving solutions to file and then comparing them using this mini-program. 
  1.4 A small matrix to be used as input, just to illustrate that the "compare_solvers" program is working. Available in /input/beamp.txt
  
  For the rest of the document, we will focus on the usage of the main program, "compare_solvers"

2. Compilation
  2.1 Dependencies:
    2.1.1 PETSc library:
    --> any recent version, personally have used 3.9.0 , also tested with PETSc 3.9.3
    --> to install PETSc library: follow instructions at https://www.mcs.anl.gov/petsc/documentation/installation.html
    --> PETSc can install all it's prerequisites by itself. In order to get access to some interesting preconditioners for testing purposes, you could install petsc with MUMPS and with HYPRE, see instructions at https://www.mcs.anl.gov/petsc/miscellaneous/external.html
    2.1.2 Intel MKL library
    --> any recent version should do, I have personally used mkl_2018.2.199
  2.2 Makefile
    --> An example makefile is provided in src/makefile, however this will most likely have to be modified
    --> the 4 variables that need to be modified from system to system inside the makefile are:
        - PETSC_ARCH , PETSC_DIR - Instructions on how to set these are available in the PETSc installation guide: https://www.mcs.anl.gov/petsc/documentation/installation.html
        - MKL_INC, MKL_SHLIB - These variables control the compiling/linking against mkl, they should be easy to set based on other makefiles of other programs that use Intel Pardiso (IE Calculix)..
   
   If dependencies were installed sucessfully and makefile was modified such that the 4 variables are correct to build use :
   
   cd src
   make release
   
   The compiled program should be the output, named "time_solvers"
  
3. Usage
  --> The behaviour of the program is fully controlled by the command line arguments. Furthermore, the PETSc library uses MPI parallelization whereas Intel Pardiso uses OpenMP parallelization. Hence the command line varies considerably depending on what exactly you want to test. Some examples are avalable in src/ (IE src/example_script_laptop example_script_slurmcluster_1 example_script_slurmcluster_2). Below a breakdown of run script is given:
  
   
   {MPI_options} ./time_solvers  ${control_string} ${inputfile_path} ${outputfile_path} ${PETSC_options}
   
   Where:
   3.1. ${control_string} is a placeholder for a string that controls the overall behaviour of the program. Possible control strings:
       3.1.1 --solve-with-pardiso
            Solves input system of linear eq using only the pardiso library, and displays (to standard output) information on the memory and time performance of the library. When this control string is used program should always be called with only 1 MPI process. The Intel Pardiso parallelization is to be controlled via an environmental variable, OMP_NUM_THREADS, that should be set before calling the program. By default solution is not printed to the output file: if this is desired, there is a line in main that can be uncommented
       3.1.2 --solve-with-PETSc
            Solves input system of linear eq using only the PETSc library, and displays (to standard output) information on the time and memory performance of the library. When this control string is used the parallelization of the PETSc library is controlled via ${MPI_options}, and the parameter OMP_NUM_THREADS has no effect. PETSc is highly configurable from the command line, in terms of what solver is used, what information is displayed etc. These configurations should be place in placeholder ${PETSC_options}, which we will talk about later. By default solution is not printed to the output file: if this is desired, there is a line in main that can be uncommented
       3.1.3 --solve-with-both
            Solves the input system with the Intel Pardiso library, followed by solving the system with the PETSc library. This option is meant to be used on input files that do not contain a solution, and provides a way to cross test the solutions of the two libraries, one against each other. Hence this option should only be used when input file does not contain a solution. In this case the parallelism of Intel Pardiso is set via OMP_NUM_THREADS, and the parallelism of PETSc is set via the MPI options: However, the MPI options should be specified so that each MPI Process runs on at least OMP_NUM_THREADS number of cores. Intel Pardiso will always run on MPI node 0 while the other nodes are waiting, followed by PETSc running on all MPI nodes specified. As in the case of --solve-with-PETSc, ${PETSC_options} can contain a myriad of options that control PETSc behaviour
       3.1.4 --create-structsym-mat and --compare-matrices
            Options which were used during development/debugging, they don't provide anything purposeful, hence should not be used
   
   3.2. ${MPI_options} is a placeholder for all settings which have to do with the MPI environment
       --> when using option --solve-with-pardiso this should be left blank or only 1 mpi process should be specified ie mpiexec.hydra -n 1
       --> for options --solve-with-PETSc or --solve-with-both multiple MPI nodes can be specified.
   3.3. ${inputfile_path} is a placeholder for a string which specifies the location of the input file, which contains
the linear system of equations to be solved. The input file must have same format as matrices provided by MTU for testing..
       --> For example ${inputfile_path} could be /home/mot_tudor/workspace/unijob/input/beamp.txt
   3.4. ${outputfile_path} is a placeholder for the path to a file that will store the solution of the solvers. Note that 
unless you modify main.cpp nothis is printed in this file, however a file should still be provided even if there is no printing involved..
   3.5. ${PETSC_options} is a placeholder for all PETSc options. For details on the PETSc options one should read the manual: http://www.mcs.anl.gov/petsc/petsc-current/docs/manual.pdf
      --> for example ${PETSC_options} could be: -mat_type seqaij  -ksp_converged_reason -pc_type ilu -pc_factor_levels 0 -pc_factor_nonzeros_along_diagonal -ksp_type gmres -pc_factor_mat_ordering_type rcm  -ksp_monitor_true_residual -ksp_gmres_restart 100 -ksp_view_pre -ksp_view -ksp_norm_type unpreconditioned -ksp_monitor_singular_value
      --> note that there are required PETSc options and optional PETSc options.. See provided script examples in src/ for more details..
  
  4. What next?
    After you have managed to succesfully run the program, with desired options (IE petsc or pardiso etc), all the data that you need for comparison of the two solvers was displayed on the standard output. Based on your environment (IE on a cluster), this standard output is a file. If it is not a file, you can easily redirect it to a file.
    The python script described in 1.2 can read these files and generate some nice graphs of the performance in time of the PETSc solver, compared to the pardiso solver (please read the separate README of the python script on how exactly to do that)


Hope this little program helps!
