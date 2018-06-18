#ifndef PETSC_SOLVER_H_INCLUDED
#define PETSC_SOLVER_H_INCLUDED

#include "petscksp.h"

//this class encapsulates the PETSc solver, and also hides the parallelism (for good or for worse)
//all parallel data structures (matrix, vectors) constructed in constructor of this class

class PETSc_solver:public lin_sys_solver{
private:
    Vec            x, b, u;      /* approx solution, RHS, exact solution */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PetscErrorCode ierr;
    PetscInt       nlocal;
    PetscInt       rstart,rend; //rowstart, rowend

    PetscErrorCode create_petsc_mat(linear_sys& input_sys);
    PetscErrorCode create_petsc_vecs(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_MPISBAIJ(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_MPIAIJ(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_SEQAIJ(linear_sys& input_sys);
    PetscErrorCode mat_preallocate_mem_SEQSBAIJ(linear_sys& input_sys);

    void investigate_paral_scaling(linear_sys& input_sys, PetscScalar time);

    //void mat_calculate_local_sizes(linear_sys & mat_dim);


public:
    virtual int solve_sys(linear_sys &sys);
    virtual void print_sol_to_file(linear_sys &sys, std::string & outputfile);
    virtual ~PETSc_solver();
    PETSc_solver(int& main_argc, char**& main_argv, linear_sys& input_sys);


    void check_petsc_solution(linear_sys &sys);
    /*provided as a hopefully better alternative to the function above, this one uses
     * PETScs VecNorm to check the relative error rather than my own calculation, which shouldn't be trusted too much*/
    void check_petsc_solution_alternative(linear_sys &sys);


    PETSc_solver() = delete;
    PETSc_solver(const linear_sys& dummy) = delete;
    const PETSc_solver& operator=(const PETSc_solver& dummy) = delete;
    PETSc_solver(const PETSc_solver&) = delete;



};

#endif // PETSC_SOLVER_H_INCLUDED
