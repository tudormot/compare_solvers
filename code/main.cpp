#include<iostream>
#include "sys_mat.h"
#include "testing.h"
#include "sys_solver.h"
#include "timer.h"

int main(int argc, char *argv[])
{

    PetscMPIInt    no_of_nodes, node_rank;
    PetscErrorCode ierr;

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &no_of_nodes);CHKERRQ(ierr);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &node_rank); CHKERRQ(ierr);


    std::string out_filename(argv[3]); //todo:could improve using move semantics but had some problems with the compilers on cluster
    std::string in_filename(argv[2]);

    linear_sys lin_sys(in_filename,no_of_nodes,node_rank);
    if(std::string(argv[1]) == "--solve-with-pardiso")
    {
        if(node_rank == 0)
        {
            if(no_of_nodes != 1)
            {
                std::cout<<"Error: pardiso is meant to be called with no_of_nodes = 1 \n";
                throw;
            }
            else
            {
                //solve with pardiso here
                pardiso_solver pard_solve(lin_sys,std::move(out_filename));
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
        PETSc_solver petsc_solve(argc_petsc,argv_petsc,lin_sys,std::move(out_filename)); //ignore the first 2 command line arguments, they are used for non PETSc stuff
        lin_sys.release_mem_mat();
        petsc_solve.solve_sys(lin_sys);
        petsc_solve.print_sol_to_file(lin_sys);
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
            pardiso_solver pard_solve(lin_sys,std::move(out_filename));
            //as this option should thoretically be used only for inputs that do not have a given solution, the pardiso solver
            //will store the solution in the lin_sys object, which is then compared agains the petsc solution in the petsc routines
            pard_solve.solve_sys(lin_sys);
    	}
    	MPI_Barrier(MPI_COMM_WORLD);
    	int argc_petsc = argc-3;
		char ** argv_petsc = argv+3;
		PETSc_solver petsc_solve(argc_petsc,argv_petsc,lin_sys,std::move(out_filename)); //ignore the first 2 command line arguments, they are used for non PETSc stuff
		lin_sys.release_mem_mat();
		petsc_solve.solve_sys(lin_sys);
		petsc_solve.print_sol_to_file(lin_sys);
    }
    else if(std::string(argv[1]) == "--create-structsym-mat")
    {
    	if(no_of_nodes != 1)
		{
			std::cout<<"Error: create-structsym-mat is meant to be called with no_of_nodes = 1 \n";
			throw;
		}
		else
		{
			test::create_structsymmat_from_symmat(lin_sys);
			//std::cout<<"DEBUG. printing original sol in "
		}
    }
    else if(std::string(argv[1]) == "--compare-matrices")
    {

    	if(no_of_nodes != 1)
		{
			std::cout<<"Error: compare-matrices is meant to be called with no_of_nodes = 1 \n";
			throw;
		}
		else
		{
			/*in this case using argv3(via out_filename) as the name of the 2nd input file*/
			linear_sys lin_sys2(out_filename,no_of_nodes,node_rank);
			test::compare_matrices(lin_sys,lin_sys2);
		}
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
