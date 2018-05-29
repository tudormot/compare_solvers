#include<iostream>
#include "sys_mat.h"
#include "testing.h"
#include "sys_solver.h"
#include "timer.h"

std::string parallel_timer::output_filename = "dummy";

int main(int argc, char *argv[])
{

    PetscMPIInt    no_of_nodes, node_rank;
    PetscErrorCode ierr;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &no_of_nodes);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &node_rank); CHKERRQ(ierr);

    if(node_rank==0)
    {
        parallel_timer::output_filename = argv[3];
    }
    std::string filename(argv[2]);
    linear_sys lin_sys(filename,no_of_nodes,node_rank);

    if(std::string(argv[1]) == "--solve-with-pardiso")
    {
        if(node_rank == 0)
        {
            if(no_of_nodes != 1)
            {
                std::cout<<"Error: pardiso is meant to be called with no_of_nodes = 1 \n";
            }
            else
            {
                //solve with pardiso here
                pardiso_solver pard_solve(lin_sys);
                lin_sys.release_mem_mat();
                pard_solve.solve_sys(lin_sys);

            }
        }
    }
    else if(std::string(argv[1]) == "--solve-with-PETSc")
    {
        //solve with PETSc
        int argc_petsc = argc-3;
        char ** argv_petsc = argv+3;
        PETSc_solver petsc_solve(argc_petsc,argv_petsc,lin_sys); //ignore the first 2 command line arguments, they are used for non PETSc stuff
        lin_sys.release_mem_mat();
        petsc_solve.solve_sys(lin_sys);
    }
    else if(std::string(argv[1]) == "--solve-with-both")
    {
    	/* solve system with pardiso first on main node, then solve with PETSc on all nodes, then compare the results of the 2*/
    	/*ideally pardiso should run on same number of OMP threads as the number of MPI processes used for PETSC.
    	 * For Example, a good combination on the mpp2 cluster is 4 mpi processes, each process having allocated 7 CPUS
    	 *  (for good socket utilization)  and each with 4 OMPthreads*/

    	/*NOTE: order of the solver matters! pardiso should be first
    	 * this is because for efficiency reasons petsc changes lin_sys vectors from 1-based indexing to 0-based indexing,
    	 * hence lin_sys can't be used to afterwards setup pardiso*/
    	if(node_rank == 0)
    	{
            pardiso_solver pard_solve(lin_sys);
            //as this option should theretically be used only for inputs that do not have a given solution, the pardiso solver
            //will store the solution in the lin_sys object, which is then compared agains the petsc solution in the petsc routines
            pard_solve.solve_sys(lin_sys);
    	}
    	MPI_Barrier(MPI_COMM_WORLD);
    	int argc_petsc = argc-3;
		char ** argv_petsc = argv+3;
		PETSc_solver petsc_solve(argc_petsc,argv_petsc,lin_sys); //ignore the first 2 command line arguments, they are used for non PETSc stuff
		lin_sys.release_mem_mat();
		petsc_solve.solve_sys(lin_sys);
    }
    else
    {
        if(node_rank == 0)
        {
            std::cout<<"unrecognised command line input\n";
        }
    }



    // Finalize the MPI environment.
    MPI_Finalize();

    return 0;
}
