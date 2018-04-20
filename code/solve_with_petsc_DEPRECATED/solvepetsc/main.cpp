#include <iostream>


#include "petscksp.h" //also include MPI api

static const char help[] = "First attempt at solving system using PETSc. Matrix is built sequencially, needs to be discussed\n\
with  MTU team whether it is feasible for the matrix to be built in parallel.";



int main(int argc, char *argv[])
{

    Vec            x, b, u;      /* approx solution, RHS, exact solution */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    PetscReal      norm;         /* norm of solution error */
    PetscErrorCode ierr;
    PetscInt       i,n = 10,col[3],its;
    PetscMPIInt    no_of_nodes, node_rank;
    PetscScalar    one = 1.0,value[3];
    PetscBool      nonzeroguess = PETSC_FALSE,changepcside = PETSC_FALSE;


    ierr = PetscInitialize(&argc, &argv,(char*)0,help); //this also initialises MPI

    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&no_of_nodes);
    CHKERRQ(ierr);
    if (no_of_nodes != 1) SETERRQ(PETSC_COMM_WORLD,1,"This is a uniprocessor example only!"); //TODO this is temp,while building it on my PC


    MPI_Comm_rank(MPI_COMM_WORLD, &node_rank); //node_rank is different for each node



    std::cout << "Hello world!ierr is "<<ierr <<"and comm rank is "<<node_rank<< std::endl;
    return 0;
}
