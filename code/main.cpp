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
