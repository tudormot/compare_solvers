#include<iostream>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include "sys_mat.h"
#include "sys_solver.h"
#include "timer.h"


int PETSc_solver::solve_sys(linear_sys& sys)
{
    //solve, with timing
    parallel_timer timer("Timing of PETSc KSPSolve routine",sys.node_rank);
    timer.start(sys.node_rank);
    KSPSolve(ksp,b,x);
    timer.stop(sys.node_rank);
    timer.display_result(sys.node_rank);
    timer.print_to_file(sys.node_rank);

    //quick test that petsc solved the system correctly
    if(sys.node_rank == 0)
    {
        std::cout<<"first elems of sol given:\n"<<sys.sol[0]<<' '<<sys.sol[1]<<' '<<sys.sol[2]<<' '<<sys.sol[3]<<'\n';
        PetscScalar sol_4elem[4];
        PetscInt sol_indices[]={0,1,2,3};
        VecGetValues(x,4,sol_indices,sol_4elem);
        std::cout<<"first elems of sol computed by PETSc:\n"<<sol_4elem[0]<<' '<<sol_4elem[1]<<' '<<sol_4elem[2]<<' '<<sol_4elem[3]<<'\n';
    }
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

    //TODO not sure if adding some elemets to the solution vector is required or not
    std::vector<PetscScalar> dummy(input_sys.mat_dim,0);
    ierr = VecSetValues(b,static_cast<PetscInt>(input_sys.mat_dim),
                 &indices[0],static_cast<PetscScalar*>(&dummy[0]),ADD_VALUES);

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
    ierr = PetscInitialize(&main_argc, &main_argv,(char*)0,help);CHKERRQ_noreturn(ierr);
    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ_noreturn(ierr);

    ierr = MatSetType(A,MATMPIAIJ);CHKERRQ_noreturn(ierr); //TODO experiment with MPISBAIJ or try the setfromarray paradigm..
    ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,input_sys.mat_dim,input_sys.mat_dim);CHKERRQ_noreturn(ierr);

    //TODO empirically tweak the parameters a bit, have some kind of preallocation...
    ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ_noreturn(ierr);
    //TODO important: this should not be required.I am not doing any preallocation which is BAD
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ_noreturn(ierr);

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
    /*ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ_noreturn(ierr);
    ierr = VecSetSizes(b,PETSC_DECIDE,input_sys.mat_dim);CHKERRQ_noreturn(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,input_sys.mat_dim);CHKERRQ_noreturn(ierr);*/



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

void PETSc_solver::view_mat()
{
    MatView(A,PETSC_VIEWER_DRAW_WORLD );
}
