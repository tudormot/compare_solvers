#ifndef PETSC_SOLVER_H_INCLUDED
#define PETSC_SOLVER_H_INCLUDED

#include "petscksp.h"

//this class encapsulates the PETSc solver, and also hides the parallelism (for good or for worse)
//all parallel data structures (matrix, vectors) constructed in constructor of this class

class PETSc_solver:public lin_sys_solver{
private:
    Vec            x, b;         /* approx solution, RHS */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PetscErrorCode ierr;		 /* used for error checking*/
    PetscInt       nlocal;       /* local number of rows, different for every mpi process*/
    PetscInt       rstart,rend;  /* rowstart, rowend, this will be different for each mpi process*/


    PetscErrorCode create_petsc_mat(linear_sys& input_sys);    //populate PETSc Mat object with matrix data stored in input_sys(data which was read from file)
    PetscErrorCode create_petsc_vecs(linear_sys& input_sys);   //preallocate memory for the PETSc Vec objects, and populate them with data stored in input _sys
    PetscErrorCode mat_preallocate_mem(linear_sys& input_sys); //preallocate memory for the PETsc Mat object

    /*different memory preallocation schemes for different PETSc internal matrix formats
     * , called from mat_preallocate_mem:*/
    PetscErrorCode mat_preallocate_mem_MPISBAIJ(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_MPIAIJ(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_SEQAIJ(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_SEQSBAIJ(linear_sys& input_sys);

    void print_ksp_info(linear_sys& input_sys, PetscScalar time);


public:

    virtual int solve_sys(linear_sys &sys); //performs the actual solving of the system, and records timing data
    virtual void print_sol_to_file(linear_sys &sys, std::string & outputfile); //prints the PETSc solution to an outputfile
    virtual ~PETSc_solver();  //deallocate all PETSc internal memory

    /*PETSc_solver(int& main_argc, char**& main_argv, linear_sys& input_sys):
     * Constructor that creates all PETSc internal data structures based on the linear system inuput stored in
     * input_sys. It also takes "commnand line" arguments using main_argc, and main_argv,
     * these parameters will influence the behaviour of the PETSc solver (which matrix type to use, which preconditioner,
     * which KSP, KSP options, logging options etc.)*/
    PETSc_solver(int& main_argc, char**& main_argv, linear_sys& input_sys);


    /* void check_petsc_solution(linear_sys &sys);
     * Function is DEPRECATED, check_petsc_solution_alternative
     * should be used instead*/
    void check_petsc_solution(linear_sys &sys);
    /*provided as a hopefully better alternative to the function above, this one uses
     * PETScs VecNorm to check the relative error rather than my own calculation*/
    void check_petsc_solution_alternative(linear_sys &sys);


    PETSc_solver() = delete;
    PETSc_solver(const linear_sys& dummy) = delete;
    const PETSc_solver& operator=(const PETSc_solver& dummy) = delete;
    PETSc_solver(const PETSc_solver&) = delete;



};

#endif // PETSC_SOLVER_H_INCLUDED
