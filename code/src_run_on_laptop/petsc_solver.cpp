#include<iostream>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include "sys_mat.h"
#include "sys_solver.h"
#include "timer.h"
#include "testing.h"


int PETSc_solver::solve_sys(linear_sys& sys)
{
    //solve, with timing
    parallel_timer timer("Timing of PETSc KSPSolve routine",sys.node_rank);
    timer.start(sys.node_rank);
    ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);
    timer.stop(sys.node_rank);
    timer.display_result(sys.node_rank);
    timer.print_to_file(sys.node_rank);

    //quick test that petsc solved the system correctly::TODO should be removed
    if(sys.node_rank == 0)
    {
        std::cout<<"first elems of sol given:\n"<<sys.sol[0]<<' '<<sys.sol[1]<<' '<<sys.sol[2]<<' '<<sys.sol[3]<<'\n';
        PetscScalar sol_4elem[4];
        PetscInt sol_indices[]={0,1,2,3};
        ierr = VecGetValues(x,4,sol_indices,sol_4elem);

        std::cout<<"first elems of sol computed by PETSc(TODO replace by better test):\n"<<sol_4elem[0]<<' '<<sol_4elem[1]<<' '<<sol_4elem[2]<<' '<<sol_4elem[3]<<'\n';
    }
    (void)(check_petsc_solution(sys));
    return true;//dummy return for the time being
}
PetscErrorCode PETSc_solver::create_petsc_mat(linear_sys& input_sys)
{
    /*NOTE:this function should always be called from only one process*/

    //first prepare the input a bit (IE input_sys).PETSc requires 0 based indexing, therefore we need to modify
    //all entries of input_sys.row_i by x = x-1;

    //TODO this could be done more efficienly in the loops that first read the file..
    for(PetscInt i = 0;i<input_sys.non_diag_no;i++)
    {
        input_sys.row_i[i]--;
    }

    for(PetscInt i = 0; i< input_sys.mat_dim-1 /*last diag elem inserted out of the loop for efficiency*/;i++)
    {
        //first add the diagonal elem
        ierr = MatSetValues(A,1,&i,1,&i,&input_sys.diag[i],ADD_VALUES);CHKERRQ(ierr);

        //now add all elems in the same row with the diagonal elem, to the right of the diagonal elem
        //and all elems in the same collumn with the diagonal elem, underneath the diagonal elem

        int row_index_start=input_sys.col_ch[i] -1;
        int row_index_end=input_sys.col_ch[i+1] -1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector

        ierr = MatSetValues(A,
                            1,&i,
                            row_index_end-row_index_start, &input_sys.row_i[row_index_start],
                            &input_sys.non_diag[row_index_start],ADD_VALUES);

        ierr = MatSetValues(A,
                            row_index_end-row_index_start, &input_sys.row_i[row_index_start],
                            1,&i,
                            &input_sys.non_diag[row_index_start],ADD_VALUES);


    }
    //add the last diag elem here
    PetscInt dummyInt = input_sys.mat_dim-1;
    ierr = MatSetValues(A,1,&dummyInt,1,&dummyInt,&input_sys.diag[dummyInt],ADD_VALUES);CHKERRQ(ierr);
    CHKERRQ(ierr);
    return ierr;

}
PetscErrorCode PETSc_solver::create_petsc_vecs(linear_sys& input_sys)
{
    //need to create an indices array containing 0,1,2..mat_dim-1 TODO is there a better api for my case?
    std::vector<PetscInt> indices(input_sys.mat_dim);
    std::iota (std::begin(indices), std::end(indices), 0);

    ierr = VecSetValues(b,static_cast<PetscInt>(input_sys.mat_dim),
                 &indices[0],static_cast<PetscScalar*>(&input_sys.rhs[0]),ADD_VALUES);

    CHKERRQ(ierr);


    return ierr;
}
PETSc_solver::~PETSc_solver()
{
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    PetscFinalize();

}
//this ugly define used as we can't have return statements in constructors
#define CHKERRQ_noreturn(ierr)          do {if (PetscUnlikely(ierr)) std::cout<<"petsc ERR! : "<<PetscError(PETSC_COMM_SELF,__LINE__,PETSC_FUNCTION_NAME,__FILE__,ierr,PETSC_ERROR_REPEAT," ")<<'\n';} while (0)

PETSc_solver::PETSc_solver(int& main_argc, char**& main_argv, linear_sys& input_sys)
{
    //todo implementing the option of using a petsc viwewr would be nice as well
    ierr = PetscInitialize(&main_argc, &main_argv,(char*)0,NULL);CHKERRQ_noreturn(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ_noreturn(ierr);

    ierr = MatSetType(A,MATMPIAIJ);CHKERRQ_noreturn(ierr); //TODO experiment with MPISBAIJ or try the setfromarray paradigm..
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,input_sys.mat_dim,input_sys.mat_dim);CHKERRQ_noreturn(ierr);

    //TODO empirically tweak the parameters a bit, have some kind of preallocation...
    ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ_noreturn(ierr);
    //TODO important: this should not be required.I am not doing any preallocation which is BAD
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ_noreturn(ierr);
    ierr = MatSetUp(A);CHKERRQ_noreturn(ierr);

    //check that the matrix is symmetric, TODO at the moment non-symmetric case not implemented
    if(input_sys.is_asymmetric == true)
    {
        PetscPrintf(PETSC_COMM_WORLD,"ERROR! asymmetric matrices are not currently supported!");
    }
    else
    {
        MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);
    }

    ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,input_sys.mat_dim,&b);CHKERRQ_noreturn(ierr);
    ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,input_sys.mat_dim,&x);CHKERRQ_noreturn(ierr);

    if(input_sys.node_rank == 0)
    {
        //sequencial part of code runnning only on main processor, whose linear_sys object contains our data of interest
        ierr = create_petsc_mat(input_sys);
        CHKERRQ_noreturn(ierr);
        ierr = create_petsc_vecs(input_sys);
        CHKERRQ_noreturn(ierr);
    }
    //now assembly , on all processors
    MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(b);
    //VecAssemblyBegin(x);
    VecAssemblyEnd(b);


    //now the solver part:
    ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
    CHKERRQ_noreturn(ierr);
    ierr = KSPSetOperators(ksp,A,A);
    CHKERRQ_noreturn(ierr);

    ierr = KSPSetFromOptions(ksp); //command line arguments required
    CHKERRQ_noreturn(ierr);

}
/*checks solution stored in the petsc parallel vector structure by performing a gather operation towards node0, where
the solution of the system, provided as input in the matrix file,is stored. A comparison between the prvided solution
and the petsc solution is made,using a "testing.h" test*/
void PETSc_solver::check_petsc_solution(linear_sys &sys)
{
    PetscInt local_size;
    PetscScalar* local_array_start;
    ierr = VecGetLocalSize(this->x,&local_size); CHKERRQ_noreturn(ierr);
    ierr = VecGetArray(this->x,&local_array_start); CHKERRQ_noreturn(ierr);

    std::cout<<"local_size_is :"<<local_size<<"\n";
    /*
    std::cout<<"debug, a simpler mpi gather exercise:\n";
    if(sys.node_rank==0)
    {
        int storage_array[sys.no_of_nodes];
        MPI_Gather(&(sys.node_rank),1,MPIU_INT,storage_array,1,MPIU_INT,0,PETSC_COMM_WORLD);

        std::cout<<"disaplying received data: "<<storage_array[0]<<'#'<<storage_array[1]<<'#'<<storage_array[2]<<'#'<<storage_array[3]<<"#\n";
    }
    else
    {
        MPI_Gather(&(sys.node_rank),1,MPIU_INT,NULL,0,MPIU_INT,0,PETSC_COMM_WORLD);
    }*/


    if(sys.node_rank == 0)
    {
        PetscScalar sol_array[sys.mat_dim];
        PetscInt    local_sizes[sys.no_of_nodes];
        PetscInt    displ[sys.no_of_nodes], rcounts[sys.no_of_nodes];

        //first must perform a gather to main node of all local sizes
        MPI_Gather(&local_size,1,MPIU_INT,local_sizes,1,MPIU_INT,0,PETSC_COMM_WORLD);

        //now initialise data to pe passed to MPIGatherv in arrays displs and rcounts
        rcounts[0]=local_sizes[0];
        displ[0]=0;
        for(int i=1;i<sys.no_of_nodes;i++)
        {
            rcounts[i]=local_sizes[i];
            displ[i]=displ[i-1]+rcounts[i];
        }

        //now get the data in the sol_array on node 0
        MPI_Gatherv(local_array_start,local_size,MPIU_SCALAR,sol_array,rcounts,displ,MPIU_SCALAR,0,PETSC_COMM_WORLD);

        //check that input file contained a solution aswell
        if(sys.sol.size()==0)
        {std::cout<<"ERR! input file did not contain a solution ! \n";}
        else
        {
            printf("testing accuracy of PETSc solution:\n");
            test::calculate_sol_tolerance(&sys.sol[0],sol_array,sys.mat_dim);
        }

    }
    else
    {
        //first must perform a gather to main node of all local sizes
        MPI_Gather(&local_size,1,MPIU_INT,NULL,1,MPIU_INT,0,PETSC_COMM_WORLD);

        //now get the data in the sol_array on node 0
        MPI_Gatherv(local_array_start,local_size,MPIU_SCALAR,NULL,NULL,NULL,MPIU_SCALAR,0,PETSC_COMM_WORLD);

    }
    VecRestoreArray(x,&local_array_start);

}

void PETSc_solver::view_mat()
{
    //MatView(A,PETSC_VIEWER_DRAW_WORLD );
}
