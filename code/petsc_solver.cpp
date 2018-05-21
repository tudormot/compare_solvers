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
	parallel_timer timer("Timing of PETSc KSPSolve routine", sys.node_rank);
	timer.start(sys.node_rank);
	ierr = KSPSolve(ksp, b, x);
	CHKERRQ(ierr);
	timer.stop(sys.node_rank);
	timer.display_result(sys.node_rank);
	timer.print_to_file(sys.node_rank);
	(void) (check_petsc_solution(sys));

	return true; //dummy return for the time being
}
PetscErrorCode PETSc_solver::create_petsc_mat(linear_sys& input_sys)
{
	/*NOTE:this function should always be called from only one process*/

	//first prepare the input a bit (IE input_sys).PETSc requires 0 based indexing, therefore we need to modify
	//all entries of input_sys.row_i by x = x-1;
	//TODO this could be done more efficiently in the loops that first read the file..
	for (PetscInt i = 0; i < input_sys.non_diag_no; i++)
	{
		input_sys.row_i[i]--;
	}

	for (PetscInt i = 0;
			i < input_sys.mat_dim - 1 /*last diag elem inserted out of the loop for efficiency*/;
			i++)
	{
		//first add the diagonal elem
		ierr = MatSetValues(A, 1, &i, 1, &i, &input_sys.diag[i], ADD_VALUES);
		CHKERRQ(ierr);

		//now add all elems in the same row with the diagonal elem, to the right of the diagonal elem
		//and all elems in the same collumn with the diagonal elem, underneath the diagonal elem

		int row_index_start = input_sys.col_ch[i] - 1;
		int row_index_end = input_sys.col_ch[i + 1] - 1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector

		ierr = MatSetValues(A, 1, &i, row_index_end - row_index_start,
				&input_sys.row_i[row_index_start],
				&input_sys.non_diag[row_index_start], ADD_VALUES);

		ierr = MatSetValues(A, row_index_end - row_index_start,
				&input_sys.row_i[row_index_start], 1, &i,
				&input_sys.non_diag[row_index_start], ADD_VALUES);

	}
	//add the last diag elem here
	PetscInt dummyInt = input_sys.mat_dim - 1;
	ierr = MatSetValues(A, 1, &dummyInt, 1, &dummyInt,
			&input_sys.diag[dummyInt], ADD_VALUES);
	CHKERRQ(ierr);CHKERRQ(ierr);
	return ierr;

}
PetscErrorCode PETSc_solver::create_petsc_vecs(linear_sys& input_sys)
{
	//need to create an indices array containing 0,1,2..mat_dim-1 TODO is there a better api for my case?
	std::vector<PetscInt> indices(input_sys.mat_dim);
	std::iota(std::begin(indices), std::end(indices), 0);

	ierr = VecSetValues(b, static_cast<PetscInt>(input_sys.mat_dim),
			&indices[0], static_cast<PetscScalar*>(&input_sys.rhs[0]),
			ADD_VALUES);

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

PETSc_solver::PETSc_solver(int& main_argc, char**& main_argv,
		linear_sys& input_sys)
{
	//todo implementing the option of using a petsc viwewr would be nice as well
	parallel_timer timer1("Timing of PETSc INIT1", input_sys.node_rank);
	timer1.start(input_sys.node_rank);
	ierr = PetscInitialize(&main_argc, &main_argv, (char*) 0, NULL);
	CHKERRCONTINUE(ierr);

	ierr = VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, input_sys.mat_dim, &b); //from what I understand this reduces to non-mpi vector if no of mpi_processes == 1
	CHKERRCONTINUE(ierr);
	//ierr = VecSetFromOptions(b);CHKERRCONTINUE(ierr);
	ierr = VecDuplicate(b, &x);
	CHKERRCONTINUE(ierr);

	ierr = VecGetOwnershipRange(x, &rstart, &rend);
	CHKERRCONTINUE(ierr);
	ierr = VecGetLocalSize(x, &nlocal);
	CHKERRCONTINUE(ierr);

	ierr = MatCreate(PETSC_COMM_WORLD, &A);
	CHKERRCONTINUE(ierr);

	ierr = MatSetFromOptions(A);
	CHKERRCONTINUE(ierr);

	ierr = MatSetSizes(A, nlocal, nlocal, input_sys.mat_dim, input_sys.mat_dim);
	CHKERRCONTINUE(ierr);

	timer1.stop(input_sys.node_rank);
	timer1.display_result(input_sys.node_rank);

	ierr = mat_preallocate_mem(input_sys);
	CHKERRCONTINUE(ierr);

	//check that the matrix is symmetric, TODO at the moment non-symmetric case not implemented
	if (input_sys.is_asymmetric == true)
	{
		PetscPrintf(PETSC_COMM_WORLD,
				"ERROR! asymmetric matrices are not currently supported!");
	}
	else
	{
		MatSetOption(A, MAT_SYMMETRIC, PETSC_TRUE);
	}

	if (input_sys.node_rank == 0)
	{
		//sequential part of code running only on main processor, whose linear_sys object contains our data of interest
		ierr = create_petsc_mat(input_sys);
		CHKERRCONTINUE(ierr);
		ierr = create_petsc_vecs(input_sys);
		CHKERRCONTINUE(ierr);
	}
	//now assembly , on all processors

	parallel_timer timer2("Timing of PETSc INIT2", input_sys.node_rank);
	timer2.start(input_sys.node_rank);
	MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
	VecAssemblyBegin(b);
	VecAssemblyBegin(x);
	VecAssemblyEnd(b);
	VecAssemblyEnd(x);

	//now the solver part:
	ierr = KSPCreate(PETSC_COMM_WORLD, &ksp);
	CHKERRCONTINUE(ierr);

	ierr = KSPSetOperators(ksp, A, A);
	CHKERRCONTINUE(ierr);

	ierr = KSPSetFromOptions(ksp); //command line arguments required
	CHKERRCONTINUE(ierr);

	ierr = KSPSetUp(ksp);
	CHKERRCONTINUE(ierr);
	timer2.stop(input_sys.node_rank);
	timer2.display_result(input_sys.node_rank);

}
/*checks solution stored in the petsc parallel vector structure by performing a gather operation towards node0, where
 the solution of the system, provided as input in the matrix file,is stored. A comparison between the provided solution
 and the petsc solution is made,using a "testing.h" test*/
void PETSc_solver::check_petsc_solution(linear_sys &sys)
{
	PetscInt local_size;
	PetscScalar* local_array_start;
	ierr = VecGetLocalSize(this->x, &local_size);
	CHKERRCONTINUE(ierr);
	ierr = VecGetArray(this->x, &local_array_start);
	CHKERRCONTINUE(ierr);

	if (sys.node_rank == 0)
	{
		PetscScalar sol_array[sys.mat_dim];
		PetscInt local_sizes[sys.no_of_nodes];
		PetscInt displ[sys.no_of_nodes], rcounts[sys.no_of_nodes];

		//first must perform a gather to main node of all local sizes
		MPI_Gather(&local_size, 1, MPIU_INT, local_sizes, 1, MPIU_INT, 0,
				PETSC_COMM_WORLD);

		//now initialise data to pe passed to MPIGatherv in arrays displs and rcounts
		rcounts[0] = local_sizes[0];
		displ[0] = 0;
		for (int i = 1; i < sys.no_of_nodes; i++)
		{
			rcounts[i] = local_sizes[i];
			displ[i] = displ[i - 1] + rcounts[i];
		}

		//now get the data in the sol_array on node 0
		MPI_Gatherv(local_array_start, local_size, MPIU_SCALAR, sol_array,
				rcounts, displ, MPIU_SCALAR, 0, PETSC_COMM_WORLD);

		//check that input file contained a solution aswell
		if (sys.sol.size() == 0)
		{
			std::cout << "ERR! input file did not contain a solution ! \n";
		}
		else
		{
			printf("testing accuracy of PETSc solution:\n");
			test::calculate_sol_tolerance(&sys.sol[0], sol_array, sys.mat_dim);
		}

	}
	else
	{
		//first must perform a gather to main node of all local sizes
		MPI_Gather(&local_size, 1, MPIU_INT, NULL, 1, MPIU_INT, 0,
				PETSC_COMM_WORLD);

		//now get the data in the sol_array on node 0
		MPI_Gatherv(local_array_start, local_size, MPIU_SCALAR, NULL, NULL,
				NULL,
				MPIU_SCALAR, 0, PETSC_COMM_WORLD);

	}
	VecRestoreArray(x, &local_array_start);

}

PetscErrorCode PETSc_solver::mat_preallocate_mem(linear_sys& sys)
{
	MatType type;
	ierr = MatGetType(A, &type);
	CHKERRQ(ierr);

	//TODO temp:
	std::cout << "matrix type is :" << type << std::endl;

	if (static_cast<std::string>(type) == "seqaij")
	{
		mat_preallocate_mem_SEQAIJ(sys);
	}
	else if (static_cast<std::string>(type) == "mpiaij")
	{
		mat_preallocate_mem_MPIAIJ(sys);
	}
	else if (static_cast<std::string>(type) == "mpisbaij")
	{
		mat_preallocate_mem_MPISBAIJ(sys);
	}
	else
	{
		std::cout << "ERROR not implemented yet\n";
	}
	return ierr;
}
PetscErrorCode PETSc_solver::mat_preallocate_mem_SEQAIJ(linear_sys& sys)
{
	//this function is only called if we are running solver on 1 MPI process
	PetscInt preallocation_data_SeqAIJ[sys.mat_dim];
	memset(preallocation_data_SeqAIJ, 0, sizeof(PetscInt) * sys.mat_dim);
	for (PetscInt i = 0; i < sys.mat_dim; i++)
	{
		preallocation_data_SeqAIJ[i] = preallocation_data_SeqAIJ[i]
				+ sys.col_ch[i + 1] - sys.col_ch[i] + 1; //+1 due to the diagonal elements
		int row_index_start = sys.col_ch[i] - 1;
		int row_index_end = sys.col_ch[i + 1] - 1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector
		for (PetscInt j = row_index_start; j < row_index_end; j++)
		{
			preallocation_data_SeqAIJ[sys.row_i[j] - 1]++;
		}
	}
	ierr = MatSeqAIJSetPreallocation(A, PETSC_DEFAULT,
			preallocation_data_SeqAIJ);
	CHKERRQ(ierr);
	return ierr;

}
PetscErrorCode PETSc_solver::mat_preallocate_mem_MPIAIJ(linear_sys& sys)
{

	PetscInt local_preallocation_recv_buffer_MPIAIJ[nlocal][2];
	if (sys.node_rank == 0)
	{
		PetscInt local_sizes[sys.no_of_nodes];
		PetscInt starting_rows[sys.no_of_nodes + 1];

		MPI_Gather(&nlocal, 1, MPIU_INT, local_sizes, 1, MPIU_INT, 0,
				PETSC_COMM_WORLD);

		PetscInt preallocation_data_MPIAIJ[sys.mat_dim][2];
		/*now we have to populate preallocation data with the data..*/
		//first set all the array to 0
		memset(preallocation_data_MPIAIJ, 0,
				sizeof(PetscInt) * 2 * sys.mat_dim);

		//now populate starting rows:
		starting_rows[0] = 0;
		for (int i = 1; i < sys.no_of_nodes + 1; i++)
		{
			starting_rows[i] = starting_rows[i - 1] + local_sizes[i - 1];
		}
		int current_node = 0;
		int first_index_not_in_diag_block = nlocal;
		for (PetscInt i = 0; i < sys.mat_dim - 1; i++)
		{
			if (i >= starting_rows[current_node + 1])
			{
				current_node++;
				first_index_not_in_diag_block += local_sizes[current_node];
			}
			preallocation_data_MPIAIJ[i][0]++; //the diag elem
			int row_index_start = sys.col_ch[i] - 1;
			int row_index_end = sys.col_ch[i + 1] - 1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector

			for (int j = row_index_start; j < row_index_end; j++)
			{
				if (sys.row_i[j] - 1 < first_index_not_in_diag_block)
				{
					preallocation_data_MPIAIJ[i][0]++;
					preallocation_data_MPIAIJ[sys.row_i[j] - 1][0]++;
				}
				else
				{
					preallocation_data_MPIAIJ[i][1]++;
					preallocation_data_MPIAIJ[sys.row_i[j] - 1][1]++;
				}
			}
		}
		preallocation_data_MPIAIJ[sys.mat_dim - 1][0]++; //last entry of diagonal
		//now send preallocation data to all nodes
		for (int i = 0; i < sys.no_of_nodes; i++)
		{
			starting_rows[i] = starting_rows[i] * 2;
			local_sizes[i] = local_sizes[i] * 2;
		}

		MPI_Scatterv(preallocation_data_MPIAIJ, local_sizes, starting_rows,
		MPIU_INT, local_preallocation_recv_buffer_MPIAIJ, nlocal * 2, MPIU_INT,
				0, PETSC_COMM_WORLD);
	}
	else
	{
		MPI_Gather(&nlocal, 1, MPIU_INT, NULL, 0, MPIU_INT, 0,
				PETSC_COMM_WORLD);
		MPI_Scatterv(NULL, NULL, NULL, MPIU_INT,
				local_preallocation_recv_buffer_MPIAIJ, nlocal * 2,
				MPIU_INT, 0, PETSC_COMM_WORLD);
	}
	//now, we have the local preallocation information on all nodes, but it is a bit jumbled up, not how PETSc wants it, need to do some array manipulation
	PetscInt d_nnz_MPIAIJ[nlocal];
	PetscInt o_nnz_MPIAIJ[nlocal];
	for (int i = 0; i < nlocal; i++)
	{
		d_nnz_MPIAIJ[i] = local_preallocation_recv_buffer_MPIAIJ[i][0];
		o_nnz_MPIAIJ[i] = local_preallocation_recv_buffer_MPIAIJ[i][1];
	}

	ierr = MatMPIAIJSetPreallocation(A, 0, d_nnz_MPIAIJ, 0, o_nnz_MPIAIJ);

	CHKERRQ(ierr);

	return ierr;
}

PetscErrorCode PETSc_solver::mat_preallocate_mem_MPISBAIJ(linear_sys& sys)
{
	PetscInt local_preallocation_recv_buffer_MPISBAIJ[nlocal][2];
	if (sys.node_rank == 0)
	{
		PetscInt local_sizes[sys.no_of_nodes];
		PetscInt starting_rows[sys.no_of_nodes + 1];
		MPI_Gather(&nlocal, 1, MPIU_INT, local_sizes, 1, MPIU_INT, 0,
				PETSC_COMM_WORLD);

		PetscInt preallocation_data_MPISBAIJ[sys.mat_dim][2];

		/*now we have to populate preallocation data with the data..*/
		//first set all the array to 0
		memset(preallocation_data_MPISBAIJ, 0,
				sizeof(PetscInt) * 2 * sys.mat_dim);

		//now populate starting rows:
		starting_rows[0] = 0;
		for (int i = 1; i < sys.no_of_nodes + 1; i++)
		{
			starting_rows[i] = starting_rows[i - 1] + local_sizes[i - 1];
		}
		int current_node = 0;
		int first_index_not_in_diag_block = nlocal;
		for (PetscInt i = 0; i < sys.mat_dim - 1; i++)
		{
			if (i >= starting_rows[current_node + 1])
			{
				current_node++;
				first_index_not_in_diag_block += local_sizes[current_node];
			}
			preallocation_data_MPISBAIJ[i][0]++; //the diag elem

			int row_index_start = sys.col_ch[i] - 1;
			int row_index_end = sys.col_ch[i + 1] - 1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector

			for (int j = row_index_start; j < row_index_end; j++)
			{
				if (sys.row_i[j] - 1 < first_index_not_in_diag_block)
				{
					preallocation_data_MPISBAIJ[i][0]++;
				}
				else
				{
					preallocation_data_MPISBAIJ[i][1]++;
				}
			}
		}
		preallocation_data_MPISBAIJ[sys.mat_dim - 1][0]++; //last entry of diagonal
		//now send preallocation data to all nodes
		for (int i = 0; i < sys.no_of_nodes; i++)
		{
			starting_rows[i] = starting_rows[i] * 2;
			local_sizes[i] = local_sizes[i] * 2;
		}

		MPI_Scatterv(preallocation_data_MPISBAIJ, local_sizes, starting_rows,
		MPIU_INT, local_preallocation_recv_buffer_MPISBAIJ, nlocal * 2,
				MPIU_INT, 0, PETSC_COMM_WORLD);
	}
	else
	{
		MPI_Gather(&nlocal, 1, MPIU_INT, NULL, 0, MPIU_INT, 0,
				PETSC_COMM_WORLD);
		MPI_Scatterv(NULL, NULL, NULL, MPIU_INT,
				local_preallocation_recv_buffer_MPISBAIJ, nlocal * 2,
				MPIU_INT, 0, PETSC_COMM_WORLD);
	}
	//now, we have the local preallocation information on all nodes, but it is a bit jumbled up, not how PETSc wants it, need to do some array manipulation
	PetscInt d_nnz_MPISBAIJ[nlocal];
	PetscInt o_nnz_MPISBAIJ[nlocal];
	for (int i = 0; i < nlocal; i++)
	{
		d_nnz_MPISBAIJ[i] = local_preallocation_recv_buffer_MPISBAIJ[i][0];
		o_nnz_MPISBAIJ[i] = local_preallocation_recv_buffer_MPISBAIJ[i][1];

	}

	ierr = MatSetBlockSizes(A,1, 1);
	CHKERRQ(ierr);
	ierr = MatMPISBAIJSetPreallocation(A, 1, 0, d_nnz_MPISBAIJ, 0,
			o_nnz_MPISBAIJ);
	CHKERRQ(ierr);


	return ierr;
}

