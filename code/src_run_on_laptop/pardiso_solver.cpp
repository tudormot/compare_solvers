#include<iostream>
#include <stdlib.h>
#include <math.h>
#include "sys_mat.h"
#include "sys_solver.h"
#include "timer.h"
#include "testing.h"



//TODO: put these in a header..this is the pardiso "header"
#ifdef __cplusplus
extern "C"{
#endif
void pardisoinit (void   *, int    *,   int *, int *, double *, int *);
void pardiso     (void   *, int    *,   int *, int *,    int *, int *,
                  double *, int    *,    int *, int *,   int *, int *,
                     int *, double *, double *, int *, double *);
void pardiso_chkmatrix  (int *, int *, double *, int *, int *, int *);
void pardiso_chkvec     (int *, int *, double *, int *);
void pardiso_printstats (int *, int *, double *, int *, int *, int *,
                           double *, int *);

#ifdef __cplusplus
}
#endif


void lin_sys_solver::CCSdiag_to_CRS(const std::vector<int> &row_i_old,const std::vector<int> &col_ch_old,const std::vector<double> &non_diag_old, const std::vector<double> &diag_old, //vectors storing input (variable name contains _old)
                              std::vector<int> &col_i_new,           std::vector<int> &row_ch_new     ,      std::vector<double> &elem_new,const bool is_asym)     //vectors storing output{
 {
    /*note: this function created three completely new vectors storing the data in the new format.. this might be a problem for input systems
      of 10GB(that we will have to deal with at some point).. currently we are still storing the data in the input format as I am not sure what format
      PETSc is using, so we might still required. TODO later in project: decide if this implementation is bad..
    */
    if(is_asym == true)
    {
        std::cout<<"ERR: asymmetric matrices not yet supported\n";throw;
    }
    else
    {

        //preallocate memory for vectors that we construct.
        col_i_new.reserve(row_i_old.size() + diag_old.size());
        row_ch_new.reserve(col_ch_old.size());
        elem_new.reserve(non_diag_old.size()+diag_old.size());
        //using the fact that the lower triagnle in CCS is given, which is equal to upper triangle in CRS
        row_ch_new.push_back(1);
        for (size_t j = 0;j<diag_old.size();j++)
        {
            elem_new.push_back(diag_old[j]);
            col_i_new.push_back(j+1);
            for(size_t i = col_ch_old[j]-1;i<col_ch_old[j+1]-1;i++)
            {
                elem_new.push_back(non_diag_old[i]);
                col_i_new.push_back(row_i_old[i]);
            }
            row_ch_new.push_back(col_i_new.size()+1);
        }

    /*note: probably the old vectors should release the memory after the new vectors are generated, as the old vectors are
    probably not needed anymore..*/
    }

}
pardiso_solver::pardiso_solver(const linear_sys &sys)
{
    //populate the vectors from data stored inside sys:
    CCSdiag_to_CRS(sys.row_i,sys.col_ch,sys.non_diag,sys.diag,this->col_i,this->row_ch,this->elem,sys.is_asymmetric);
}

int pardiso_solver::solve_sys(linear_sys &sys)
{
    int      nnz = static_cast<int>(elem.size());   //no of non zero elem

    int mtype = 99;//illegal type
    if(sys.is_asymmetric == true)
    {
        std::cout<<"ERR: asymmetric matrices not supported yet\n";
        throw;
        //TODO
    }
    else
    {
        mtype = -2;        /* Real symmetric matrix */
    }

    int      nrhs = 1;          /* Number of right hand sides. */

    /* Internal solver memory pointer pt,                  */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
    /* or void *pt[64] should be OK on both architectures  */
    void    *pt[64];

    /* Pardiso control parameters. */
    int      iparm[64];
    double   dparm[64];
    int      maxfct, mnum, phase, error, msglvl, solver;

    /* Number of processors. */
    int      num_procs;

    /* Auxiliary variables. */
    char    *var;
    int      i, k;

    double   ddum;              /* Double dummy */
    int      idum;              /* Integer dummy. */

    int n = static_cast<int>(sys.mat_dim);//TODO not the case anymore, change this to not duplicate data

/* -------------------------------------------------------------------- */
/* ..  Setup Pardiso control parameters.                                */
/* -------------------------------------------------------------------- */

    error = 0;
    solver=0;/* use sparse direct solver */
    pardisoinit (pt,  &mtype, &solver, iparm, dparm, &error);

    if (error != 0)
    {
        if (error == -10 )
           printf("No license file found \n");
        if (error == -11 )
           printf("License is expired \n");
        if (error == -12 )
           printf("Wrong username or hostname \n");
        throw;
    }
    else
        printf("[PARDISO]: License check was successful ... \n");

    /* Numbers of processors, value of OMP_NUM_THREADS */
    var = getenv("OMP_NUM_THREADS");
    if(var != NULL)
        sscanf( var, "%d", &num_procs );
    else {
        printf("Set environment OMP_NUM_THREADS to 1");
        exit(1);
    }
    iparm[2]  = num_procs;

    maxfct = 1;		/* Maximum number of numerical factorizations.  */
    mnum   = 1;         /* Which factorization to use. */

    msglvl = 1;         /* Print statistical information  */
    error  = 0;         /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
/*     notation.                                                        */
/* -------------------------------------------------------------------- */
//not required as our input is already in fortran 1-based notation

#ifndef INTEL_PARDISO
/* -------------------------------------------------------------------- */
/*  .. pardiso_chk_matrix(...)                                          */
/*     Checks the consistency of the given matrix.                      */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    pardiso_chkmatrix  (&mtype, &n, &elem[0], &row_ch[0], &col_i[0], &error);
    if (error != 0) {
        printf("\nERROR in consistency of matrix: %d", error);
        exit(1);
    }
    else
    {
        std::cout<<"matrix check appears to be successful..\n";
    }

/* -------------------------------------------------------------------- */
/* ..  pardiso_chkvec(...)                                              */
/*     Checks the given vectors for infinite and NaN values             */
/*     Input parameters (see PARDISO user manual for a description):    */
/*     Use this functionality only for debugging purposes               */
/* -------------------------------------------------------------------- */

    pardiso_chkvec (&n, &nrhs, const_cast<double*>(&sys.rhs[0]), &error);
    if (error != 0) {
        printf("\nERROR  in right hand side: %d", error);
        exit(1);
    }

/* -------------------------------------------------------------------- */
/* .. pardiso_printstats(...)                                           */
/*    prints information on the matrix to STDOUT.                       */
/*    Use this functionality only for debugging purposes                */
/* -------------------------------------------------------------------- */

    pardiso_printstats (&mtype, &n, &elem[0], &row_ch[0], &col_i[0], &nrhs, const_cast<double*>(&sys.rhs[0]), &error);
    if (error != 0) {
        printf("\nERROR in printstats %d", error);
        exit(1);
    }

#endif // INTEL_PARDISO

/* -------------------------------------------------------------------- */
/* ..  Now perform Analysis, numerical factorization, solve, iterative refinement, all in one pardiso call
       NOTE: Calculix does not do this , it rather does this via 2 calls: TODO later: test time spent doing that*/
/* -------------------------------------------------------------------- */
    phase = 13;

    this->sol.resize(sys.mat_dim); //sol will store the solution

    //solve, with timing
    parallel_timer timer("Timing of Pardiso's solving routine",sys.node_rank);
    timer.start(sys.node_rank);

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
	     &n, &(this->elem[0]), &(this->row_ch[0]), &(this->col_i[0]), &idum, &nrhs,
             iparm, &msglvl, const_cast<double*>(&sys.rhs[0]), &(this->sol)[0], &error, dparm);

    timer.stop(sys.node_rank);
    timer.display_result(sys.node_rank);
    timer.print_to_file(sys.node_rank);


    if (error != 0) {
        printf("\nERROR during pardiso solving: %d", error);
        exit(1);
    }
    printf("\nSolving completd... ");





/* -------------------------------------------------------------------- */
/* ..  Termination and release of memory.                               */
/* -------------------------------------------------------------------- */
    phase = -1;                 /* Release internal memory. */

    pardiso (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, &row_ch[0], &col_i[0], &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error,  dparm);

    if(sys.sol.size() != 0)
    {
        //this means that the sys object  already stores a solution, this is either due to another solver object solving the system beforehand, or the input file containing the solution as well for testing purposes
        //in this case , nothing to do, the solution just generated by the pardiso object will be stored internally in the pardiso_solver object
        std::cout<<"It appears that the system already contains a solution. Solution stored inside pardiso_solver object \n";
        printf("Relative error of solution vs solution given:\n");
        test::calculate_sol_tolerance(&sys.sol[0],&(this->sol[0]),sys.mat_dim);
    }
    else
    {
        //store the solution inside the sys object using move semantics
        std::cout<<"Storing solution in sys obj.\n";
        std::cout<<"Can't calculate relative errror in solution as err is not given\n";
        sys.sol = std::move(this->sol);
    }

    //TODO dummy return for the time being, to be used for error handling
    return 0;

}

pardiso_solver::~pardiso_solver()
{
    //std::cout<<"dummy\n";
}
lin_sys_solver::~lin_sys_solver()
{

}

