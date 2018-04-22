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

    ierr = VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,input_sys.mat_dim,&b);CHKERRQ_noreturn(ierr);
    //ierr = VecSetFromOptions(b);CHKERRQ_noreturn(ierr);
    ierr = VecDuplicate(b,&x);CHKERRQ_noreturn(ierr);

    ierr = VecGetOwnershipRange(x,&rstart,&rend);CHKERRCONTINUE(ierr);
    ierr = VecGetLocalSize(x,&nlocal);CHKERRCONTINUE(ierr);


    ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ_noreturn(ierr);
    ierr = MatSetType(A,MATMPIAIJ);CHKERRQ_noreturn(ierr); //TODO experiment with MPISBAIJ or try the setfromarray paradigm..

    ierr = MatSetSizes(A,nlocal,nlocal,input_sys.mat_dim,input_sys.mat_dim);CHKERRQ_noreturn(ierr);

    #ifdef old
    //TODO empirically tweak the parameters a bit, have some kind of preallocation...
    ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL);CHKERRQ_noreturn(ierr);
    //TODO important: this should not be required.I am not doing any preallocation which is BAD
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);CHKERRQ_noreturn(ierr);
    ierr = MatSetUp(A);CHKERRQ_noreturn(ierr);
    #endif
    ierr = mat_preallocate_mem(input_sys);CHKERRQ_noreturn(ierr);
    //check that the matrix is symmetric, TODO at the moment non-symmetric case not implemented
    if(input_sys.is_asymmetric == true)
    {
        PetscPrintf(PETSC_COMM_WORLD,"ERROR! asymmetric matrices are not currently supported!");
    }
    else
    {
        MatSetOption(A,MAT_SYMMETRIC,PETSC_TRUE);
    }



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

PetscErrorCode PETSc_solver::mat_preallocate_mem(linear_sys& sys)
{
    PetscInt local_sz[2]; //local_sz[0] is number of rows local_sz[1] is number of cols

    ierr = MatGetLocalSize(A,&local_sz[0],&local_sz[1]); CHKERRQ(ierr);
    std::cout<<"debug, local sizes are "<<local_sz[0]<<'#'<<local_sz[1]<<"in node rank "<<sys.node_rank<<'\n';

    PetscInt local_preallocation_recv_buffer[local_sz[0]][2];
    if(sys.node_rank == 0)
    {
        PetscInt local_sizes[sys.no_of_nodes][2];
        PetscInt starting_rows[sys.no_of_nodes+1];

        MPI_Gather(local_sz,2,MPIU_INT,local_sizes,2,MPIU_INT,0,PETSC_COMM_WORLD);

        PetscInt preallocation_data[sys.no_of_nodes * local_sz[0]][2]; //last few entries might not get populated, but that is ok
        //TODO temp sanity check that local number of colums is equal:
        for(int i=0;i<sys.no_of_nodes-1;i++)
        {
            if(local_sizes[i][1] != local_sizes[i+1][1])
            {
                std::cout<<"ERR IN PREALLOCATE MAT: local number of columns not equal!\n";
            }
        }//TODO to be removed..
#if 1
        /*now we have to populate preallocation data with the data..*/
        //first set all the array to 0
        memset(preallocation_data,0,sizeof(PetscInt)*2*sys.no_of_nodes*local_sz[0]);
        //now populate starting rows:
        starting_rows[0]=0;
        for (int i =1 ;i<sys.no_of_nodes+1;i++)
        {
            starting_rows[i]=starting_rows[i-1]+local_sizes[i-1][0];
        }
        //now preallocation data..
        int current_node=0;
        int first_index_not_in_diag_block = local_sz[0];
        for(PetscInt i = 0; i< sys.mat_dim-1;i++)
        {
            if(i>=starting_rows[current_node+1])
            {
                current_node++;
                first_index_not_in_diag_block+=local_sz[0];
            }
            preallocation_data[i][0]++; //the diag elem
            int row_index_start=sys.col_ch[i] -1;
            int row_index_end=sys.col_ch[i+1] -1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector

            for(int j = row_index_start;j<row_index_end;j++)
            {
                if(sys.row_i[j]<first_index_not_in_diag_block)
                {
                    preallocation_data[i][0]++;
                    preallocation_data[sys.row_i[j]-1][0]++;
                }
                else
                {
                    preallocation_data[i][1]++;
                    preallocation_data[sys.row_i[j]-1][1]++;
                }
            }
        }
        preallocation_data[sys.mat_dim-1][0]++; //last entry of diagonal

        //now send preallocation data to all nodes
        MPI_Scatter(preallocation_data,local_sz[0]*2,MPIU_INT,local_preallocation_recv_buffer,local_sz[0]*2,MPIU_INT,0,PETSC_COMM_WORLD);
#endif
    }
    else
    {
        MPI_Gather(local_sz,2,MPIU_INT,NULL,0,MPIU_INT,0,PETSC_COMM_WORLD);
        #if 1
        MPI_Scatter(NULL,0,MPIU_INT,local_preallocation_recv_buffer,local_sz[0]*2,MPIU_INT,0,PETSC_COMM_WORLD);
        #endif
    }
    //now, we have the local preallocation information on all nodes, but it is a bit jumbled up, not how PETSc wants it, need to do some array manipulation
    #if 1
    PetscInt d_nnz[local_sz[0]];
    PetscInt o_nnz[local_sz[0]];
    for(int i=0;i<local_sz[0];i++)
    {
        d_nnz[i]=local_preallocation_recv_buffer[i][0];
        o_nnz[i]=local_preallocation_recv_buffer[i][1];
    }
    ierr = MatMPIAIJSetPreallocation(A,0,d_nnz,0,o_nnz); CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);//todo TEMP
    #else
    ierr = MatMPIAIJSetPreallocation(A,5,NULL,5,NULL); CHKERRQ(ierr);
    ierr = MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE); CHKERRQ(ierr);

    std::cout<<"debug, local sizes are "<<local_sz[0]<<'#'<<local_sz[1]<<"in node rank "<<sys.node_rank<<'\n';

    #endif // debug

    return 0;
}
#if old
void PETSc_solver::mat_calculate_local_sizes(linear_sys& sys)
{
    if(sys.mat_dim%sys.no_of_nodes==0)
    {
        this->nrows_loc = sys.mat_dim/sys.no_of_nodes;
    }
    else
    {
        if(sys.node_rank == sys.no_of_nodes-1)
        {
            this->nrows_loc=sys.mat_dim%sys.no_of_nodes;
        }
        else
        {
            this->nrows_loc = sys.mat_dim/sys.no_of_nodes;
        }
    }
    this->ncols_loc=sys.mat_dim;
}
#endif
