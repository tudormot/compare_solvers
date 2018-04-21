#ifndef PETSC_SOLVER_H_INCLUDED
#define PETSC_SOLVER_H_INCLUDED

#include "petscksp.h"

//this class encapsulates the PETSc solver, and also hides the parallelism (for good or for worse)
//all parallel data structures (matrix, vectors) constructed in constructor of this class

class PETSc_solver:public lin_sys_solver{
private:
    //TODO:cleanup these fields, were copied and some might not be needed
    Vec            x, b, u;      /* approx solution, RHS, exact solution */
    Mat            A;            /* linear system matrix */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    PetscReal      norm;         /* norm of solution error */
    PetscErrorCode ierr;
    PetscInt       i,n = 10,col[3],its;
    PetscScalar    one = 1.0,value[3];
    PetscBool      nonzeroguess = PETSC_FALSE,changepcside = PETSC_FALSE;

    PetscErrorCode create_petsc_mat(linear_sys& input_sys);
    PetscErrorCode create_petsc_vecs(linear_sys& input_sys);



public:
    virtual int solve_sys(linear_sys &sys);
    virtual ~PETSc_solver();
    PETSc_solver(int& main_argc, char**& main_argv, linear_sys& input_sys);
    void view_mat();
    void check_petsc_solution(linear_sys &sys);

    PETSc_solver() = delete;
    PETSc_solver(const linear_sys& dummy) = delete;
    const PETSc_solver& operator=(const PETSc_solver& dummy) = delete;
    PETSc_solver(const PETSc_solver&) = delete;
};

#endif // PETSC_SOLVER_H_INCLUDED
