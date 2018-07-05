#include<iostream>
#include <algorithm>    // std::max
#include <stdlib.h>
#include <math.h>
#include "sys_mat.h"
#include "sys_solver.h"
#include "timer.h"
#include "testing.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#include "mkl_service.h"

#define CTRL_PARAM_SZ (128)

//TODO: lin_sys_solver functions should be part od a separate sys_mat.c file
void lin_sys_solver::CLCXformats_to_CRS(const std::vector<int> &row_i_old,const std::vector<int> &col_ch_old,const std::vector<double> &non_diag_old, const std::vector<double> &diag_old, //vectors storing input (variable name contains _old)
                                    std::vector<int> &col_i_new,           std::vector<int> &row_ch_new,      std::vector<double> &elem_new,const bool is_asym)          //vectors storing output{
{
    /*note: this function created three completely new vectors storing the data in the new format.. this might be a problem for input systems
      of 10GB(that we will have to deal with at some point).. currently we are still storing the data in the input format as I am not sure what format
      PETSc is using, so we might still required. TODO later in project: decide if this implementation is bad..
    */

	//temp, quick test:
	//std::cout<<"In functionc CRSdiag_to CRS. Size of input vectors\n";
	//std::cout<<row_i_old.size()<<' '<<col_ch_old.size()<<' '<<non_diag_old.size()<<' '<<diag_old.size()<<std::endl;

    if(is_asym == true)
    {
        col_i_new.resize(2 * row_i_old.size() + diag_old.size());
		row_ch_new.reserve(col_ch_old.size());
		elem_new.resize(non_diag_old.size()+diag_old.size());
		/*we assume that non_diag_old contains lower triangular stored column by colum, followed by uper triangular stored
		 * row by row*/


		/*we need some helper vectors to perform the conversion. One is a vector that initially: specifies at which indexes of elem_new
		 * the col/ row changes. To generate this vect we need another vector which speicifies how many elements are present in each row*/
		//actually these vectors can be merged in one, for efficiency reasons
		std::vector<int> helper_vect(diag_old.size()+1,0);

		//now go through input vectors once and calculate how many elements are present per matrix row, store in helper_vect
		//keep in mind that input is in 1-based indexing
		//after this for loop helper_vect will store 0 in helper_vect[0], and number of elems on row n-1 in entry n
		for (int i = 1; i < diag_old.size()+1; i++)
		{
			helper_vect[i] =helper_vect[i] + col_ch_old[i] - col_ch_old[i-1] + 1; //+1 due to the diagonal elements
			int row_index_start = col_ch_old[i-1] - 1;
			int row_index_end = col_ch_old[i] - 1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector
			for (int j = row_index_start; j < row_index_end; j++)
			{
				helper_vect[row_i_old[j]]++;
			}
		}
		//now calculate at which indices of elem_new will the row change
		for (int i =1;i<diag_old.size()+1;i++)
		{
			helper_vect[i]=helper_vect[i]+helper_vect[i-1];
		}
		//now with the help of helper vect generate row_ch new
		for(int i =0;i<diag_old.size()+1;i++)
		{
			row_ch_new.push_back(helper_vect[i]+1); // plus one because of row based indexing
		}

		//now generate col_i
		for(int i = 0; i< diag_old.size()-1; i++)
		{
			//store index of diagonal and diag elem in right location
			col_i_new[helper_vect[i]]=i+1;
			elem_new[helper_vect[i]] = diag_old[i];
			helper_vect[i]++;

			int row_index_start = col_ch_old[i] - 1;
			int row_index_end = col_ch_old[i+1] - 1; //these two variables will store the start and end indexes of the data we want to extract from linear_sys.elem vector
			for (int j = row_index_start; j < row_index_end; j++)
			{
				//now store indices of elems in line i+1, on the right of the diagonal, plus elems
				col_i_new[helper_vect[i]] = row_i_old[j];
				elem_new[helper_vect[i]]  = non_diag_old [j+non_diag_old.size()/2]; //TODO: could be written better by modifying fun param
				helper_vect[i]++;
				//now store indices of elems in col i+1, underneath the diagonal, plus elems
				col_i_new[helper_vect[row_i_old[j]-1]] = i+1; //+1 due to one based indexing
				elem_new[helper_vect[row_i_old[j]-1]]  = non_diag_old [j];
				helper_vect[row_i_old[j]-1]++;
			}
		}
		//add last diagonal entry:
		col_i_new[helper_vect[diag_old.size()-1]] = diag_old.size();
		elem_new[helper_vect[diag_old.size()-1]] = diag_old.back();


    }
    else
    {
        //preallocate memory for vectors that we construct.
        col_i_new.reserve(row_i_old.size() + diag_old.size());
        row_ch_new.reserve(col_ch_old.size());
        elem_new.reserve(non_diag_old.size()+diag_old.size());
        //using the fact that the lower triagnle in CCS is given, which is equal to upper triangle in CRS
        row_ch_new.push_back(1);
        for (size_t j = 0; j<diag_old.size(); j++)
        {
            elem_new.push_back(diag_old[j]);
            col_i_new.push_back(j+1);
            for(size_t i = col_ch_old[j]-1; i<col_ch_old[j+1]-1; i++)
            {
                elem_new.push_back(non_diag_old[i]);
                col_i_new.push_back(row_i_old[i]);
            }
            row_ch_new.push_back(col_i_new.size()+1);
        }
    }

}
void pardiso_solver::print_sol_to_file(linear_sys &sys, std::string & outputfile)
{
	std::cout<<"DUMMy not implemented";
}
pardiso_solver::pardiso_solver(const linear_sys &sys)
{
    //populate the vectors from data stored inside sys:
	CLCXformats_to_CRS(sys.row_i,sys.col_ch,sys.non_diag,sys.diag,this->col_i,this->row_ch,this->elem,sys.is_asymmetric);
}

int pardiso_solver::solve_sys(linear_sys &sys)
{
    /* Matrix data. */
    MKL_INT n = (MKL_INT)sys.mat_dim;
    MKL_INT *ia=&(this->row_ch)[0];
    MKL_INT *ja=&(this->col_i)[0];
    double *a = &(this->elem)[0];

    MKL_INT mtype=-99;//illegal magic number
    if(sys.is_asymmetric == true)
    {
        mtype =  1;       /*real structurally symmetric matrix*/
    }
    else
    {
        mtype = -2;       /* Real symmetric matrix */
    }



    /* RHS and solution vectors. */
    //rhs already stored in sys object
    this->sol.resize(sys.mat_dim); //sol will store the solution

    MKL_INT nrhs = 1;     /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i;
    double ddum;          /* Double dummy */
    MKL_INT idum;         /* Integer dummy. */
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;         /* No solver default */
    iparm[1] = 2;         /* Fill-in reordering from METIS */
    iparm[3] = 0;         /* No iterative-direct algorithm */
    iparm[4] = 0;         /* No user fill-in reducing permutation */
    iparm[5] = 0;         /* Write solution into x */
    iparm[6] = 0;         /* Not in use */
    iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    iparm[8] = 0;         /* Not in use */
    iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;        /* Not in use */
    iparm[12] = 0;        /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1 in case of inappropriate accuracy */
    iparm[13] = 0;        /* Output: Number of perturbed pivots */
    iparm[14] = 0;        /* Not in use */
    iparm[15] = 0;        /* Not in use */
    iparm[16] = 0;        /* Not in use */
    iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1;       /* Output: Mflops for LU factorization */
    iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    maxfct = 1;           /* Maximum number of numerical factorizations. */
    mnum = 1;         /* Which factorization to use. */
    msglvl = 1;           /* Print statistical information in file */
    error = 0;            /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for ( i = 0; i < 64; i++ )
    {
        pt[i] = 0;
    }
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during symbolic factorization: %d", error);
        exit (1);
    }
    printf ("\nReordering completed ... ");
    printf ("Peak memory on symbolic factorization: %d\n",iparm[14]);
    printf ("Permanent memory on symbolic factorization: %d\n",iparm[15]);
    printf ("Peak memory on numerical factorization and solution: %d\n",iparm[16]);
    printf ("Peak total memory used by pardiso (maximum of above 3 according to intel manual): %d ",std::max(iparm[14],std::max(iparm[15],iparm[16])));
    printf ("\nNumber of nonzeros in factors = %d", iparm[17]);
    printf ("\nNumber of factorization MFLOPS = %d", iparm[18]);

    /* -------------------------------------------------------------------- */
    /* .. Numerical factorization. */
    /* -------------------------------------------------------------------- */
    phase = 22;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
    if ( error != 0 )
    {
        printf ("\nERROR during numerical factorization: %d", error);
        exit (2);
    }
    printf ("\nFactorization completed ... ");
    /* -------------------------------------------------------------------- */
    /* .. Back substitution and iterative refinement. */
    /* -------------------------------------------------------------------- */
    phase = 33;
    iparm[7] = 2;         /* Max numbers of iterative refinement steps. */

    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &(sys.rhs)[0], &(this->sol)[0], &error);
    if ( error != 0 )
    {
        printf ("\nERROR during solution: %d", error);
        exit (3);
    }
    printf ("\nSolve completed ... ");
    if(sys.sol.size() != 0)
    {
	std::cout<<"size of already stored sol vector is :"<<sys.sol.size()<<'\n';
	std::cout<<"mat_dim is :"<<sys.mat_dim<<'\n';        
//this means that the sys object  already stores a solution, this is either due to another solver object solving the system beforehand, or the input file containing the solution as well for testing purposes
        //in this case , nothing to do, the solution just generated by the pardiso object will be stored internally in the pardiso_solver object
        std::cout<<"It appears that the system already contains a solution. \n";
        printf("Relative error of solution vs solution given:\n");
        test::calculate_sol_tolerance(&sys.sol[0],&(this->sol[0]),sys.mat_dim);
        this->sol.resize(0);
        this->sol.shrink_to_fit();
    }
    else
    {
        //store the solution inside the sys object using move semantics
        std::cout<<"Solution not provided in input file, storing pardiso solution as system solution\n";
        sys.sol = std::move(this->sol);
    }
    /* -------------------------------------------------------------------- */
    /* .. Termination and release of memory. */
    /* -------------------------------------------------------------------- */
    phase = -1;           /* Release internal memory. */
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
             &n, &ddum, ia, ja, &idum, &nrhs,
             iparm, &msglvl, &ddum, &ddum, &error);

    /*also freeinternal representation of matrix, as it is not needed anymore*/
    this->col_i.resize(0);
    this->col_i.shrink_to_fit();
    this->elem.resize(0);
	this->elem.shrink_to_fit();
	this->row_ch.resize(0);
	this->row_ch.shrink_to_fit();
    return 0;
}

pardiso_solver::~pardiso_solver()
{
    //std::cout<<"dummy\n";
}
lin_sys_solver::~lin_sys_solver()
{

}

int pardiso_solver::solve_sys_gmres(linear_sys &sys)
{
    #ifdef not_implemented
    if(sys.is_asymmetric == true)
    {
        std::cout<<"only symmetric matrices supported for the time being"
        goto FAILED;
    }
    /*---------------------------------------------------------------------------
    * Define arrays for the upper triangle of the coefficient matrix
    * Compressed sparse row storage is used for sparse representation
    *---------------------------------------------------------------------------*/
    MKL_INT* ia = static_cast<MKL_INT*>(&(row_ch)[0]);
    MKL_INT* ja = static_cast<MKL_INT*>(&(col_i)[0]);
    double* A = static_cast<double*>(&(elem)[0]);
    //int N = sys.mat_dim; //replacement for a previously used macro

    /*---------------------------------------------------------------------------
    * Allocate storage for the ?par parameters and the solution/rhs/residual vectors
    *---------------------------------------------------------------------------*/
    MKL_INT ipar[CTRL_PARAM_SZ];
    double dpar[CTRL_PARAM_SZ], tmp[sys.mat_dim * (2 * sys.mat_dim + 1) + (sys.mat_dim * (sys.mat_dim + 9)) / 2 + 1];
    //double expected_solution[N] = { -1.0, 1.0, 0.0, 1.0, -1.0 }; //prob not needed
    double rhs[sys.mat_dim]; //will copy rhs stored in sys object, as this solver destroyes the rhs array it uses
    //b[N]; //probably not needed
    double computed_solution[sys.mat_dim];
    double residual[sys.mat_dim];
    /*---------------------------------------------------------------------------
    * Some additional variables to use with the RCI (P)FGMRES solver
    *---------------------------------------------------------------------------*/
    MKL_INT itercount, expected_itercount = 4;
    MKL_INT RCI_request, i, ivar;
    double dvar;
    char cvar;

    printf ("--------------------------------------------------------\n");
    printf ("The FULLY ADVANCED example of usage of RCI FGMRES solver\n");
    printf ("   to solve a non-symmetric indefinite non-degenerate\n");
    printf ("          algebraic system of linear equations\n");
    printf ("--------------------------------------------------------\n\n");
    /*---------------------------------------------------------------------------
    * Initialize variables and the right hand side through matrix-vector product
    *---------------------------------------------------------------------------*/
    ivar = sys.mat_dim;
    cvar = 'N';
    std::cout<<"DEBUG, the N macro problem:"
    if(cvar==33)
    {
        std::cout<<"seems that macros get expanded inside char literals aswell.\n";
    }
    else
    {
        std::cout<<"seems that macros do not get expanded inside char literals\n";
    }

    //mkl_dcsrgemv (&cvar, &ivar, A, ia, ja, expected_solution, rhs); //not needed as we already have a rhs
    /*---------------------------------------------------------------------------
    * copy the sys solution in array rhs
    *---------------------------------------------------------------------------*/
    i = 1;
    dcopy (&ivar, &(sys.rhs)[0], &i, rhs, &i);
    /*---------------------------------------------------------------------------
    * Initialize the initial guess
    *---------------------------------------------------------------------------*/
    for (i = 0; i < N; i++)
    {
        computed_solution[i] = 0.0;
    }
    /*---------------------------------------------------------------------------
    * Initialize the solver
    *---------------------------------------------------------------------------*/
    dfgmres_init (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request != 0)
        goto FAILED;
    /*---------------------------------------------------------------------------
    * Set the desired parameters:
    * do the restart after 2 iterations
    * LOGICAL parameters:
    * do not do the stopping test for the maximal number of iterations
    * do the Preconditioned iterations of FGMRES method
    * DOUBLE PRECISION parameters
    * set the relative tolerance to 1.0D-3 instead of default value 1.0D-6

    TODO these should be modified
    *---------------------------------------------------------------------------*/
    ipar[14] = 2;
    ipar[7] = 0;
    ipar[10] = 1;
    dpar[0] = 1.0E-3;
    /*---------------------------------------------------------------------------
    * Check the correctness and consistency of the newly set parameters
    *---------------------------------------------------------------------------*/
    dfgmres_check (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar,
                   tmp);
    if (RCI_request != 0)
        goto FAILED;
    /*---------------------------------------------------------------------------
    * Print the info about the RCI FGMRES method
    *---------------------------------------------------------------------------*/
    printf ("Some info about the current run of RCI FGMRES method:\n\n");
    if (ipar[7])
    {
        printf ("As ipar[7]=%d, the automatic test for the maximal number of ", ipar[7]);
        printf ("iterations will be\nperformed\n");
    }
    else
    {
        printf ("As ipar[7]=%d, the automatic test for the maximal number of ", ipar[7]);
        printf ("iterations will be\nskipped\n");
    }
    printf ("+++\n");
    if (ipar[8])
    {
        printf ("As ipar[8]=%d, the automatic residual test will be performed\n", ipar[8]);
    }
    else
    {
        printf ("As ipar[8]=%d, the automatic residual test will be skipped\n", ipar[8]);
    }
    printf ("+++\n");
    if (ipar[9])
    {
        printf ("As ipar[9]=%d, the user-defined stopping test will be ", ipar[9]);
        printf ("requested via\nRCI_request=2\n");
    }
    else
    {
        printf ("As ipar[9]=%d, the user-defined stopping test will not be ", ipar[9]);
        printf ("requested, thus,\nRCI_request will not take the value 2\n");
    }
    printf ("+++\n");
    if (ipar[10])
    {
        printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will be ", ipar[10]);
        printf ("performed, thus,\nthe preconditioner action will be requested via ");
        printf ("RCI_request=3\n");
    }
    else
    {
        printf ("As ipar[10]=%d, the Preconditioned FGMRES iterations will not ", ipar[10]);
        printf ("be performed,\nthus, RCI_request will not take the value 3\n");
    }
    printf ("+++\n");
    if (ipar[11])
    {
        printf ("As ipar[11]=%d, the automatic test for the norm of the next ", ipar[11]);
        printf ("generated vector is\nnot equal to zero up to rounding and ");
        printf ("computational errors will be performed,\nthus, RCI_request will not ");
        printf ("take the value 4\n");
    }
    else
    {
        printf ("As ipar[11]=%d, the automatic test for the norm of the next ", ipar[11]);
        printf ("generated vector is\nnot equal to zero up to rounding and ");
        printf ("computational errors will be skipped,\nthus, the user-defined test ");
        printf ("will be requested via RCI_request=4\n");
    }
    printf ("+++\n\n");
    /*---------------------------------------------------------------------------
    * Compute the solution by RCI (P)FGMRES solver with preconditioning
    * Reverse Communication starts here
    *---------------------------------------------------------------------------*/
ONE:
    dfgmres (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
    /*---------------------------------------------------------------------------
    * If RCI_request=0, then the solution was found with the required precision
    *---------------------------------------------------------------------------*/
    if (RCI_request == 0)
        goto COMPLETE;
    /*---------------------------------------------------------------------------
    * If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
    * and put the result in vector tmp[ipar[22]-1]
    *---------------------------------------------------------------------------
    * NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
    * therefore, in C code it is required to subtract 1 from them to get C style
    * addresses
    *---------------------------------------------------------------------------*/
    if (RCI_request == 1)
    {
        mkl_dcsrgemv (&cvar, &ivar, A, ia, ja, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
        goto ONE;
    }
    /*---------------------------------------------------------------------------
    * If RCI_request=2, then do the user-defined stopping test
    * The residual stopping test for the computed solution is performed here
    *---------------------------------------------------------------------------
    * NOTE: from this point vector b[N] is no longer containing the right-hand
    * side of the problem! It contains the current FGMRES approximation to the
    * solution. If you need to keep the right-hand side, save it in some other
    * vector before the call to dfgmres routine. Here we saved it in vector
    * rhs[N]. The vector b is used instead of rhs to preserve the
    * original right-hand side of the problem and guarantee the proper
    * restart of FGMRES method. Vector b will be altered when computing the
    * residual stopping criterion!
    *---------------------------------------------------------------------------*/
    if (RCI_request == 2)
    {
        /* Request to the dfgmres_get routine to put the solution into b[N] via ipar[12]
           --------------------------------------------------------------------------------
           WARNING: beware that the call to dfgmres_get routine with ipar[12]=0 at this
           stage may destroy the convergence of the FGMRES method, therefore, only
           advanced users should exploit this option with care */
        ipar[12] = 1;
        /* Get the current FGMRES solution in the vector b[N] */
        dfgmres_get (&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
        /* Compute the current true residual via MKL (Sparse) BLAS routines */
        mkl_dcsrgemv (&cvar, &ivar, A, ia, ja, b, residual);
        dvar = -1.0E0;
        i = 1;
        daxpy (&ivar, &dvar, rhs, &i, residual, &i);
        dvar = dnrm2 (&ivar, residual, &i);
        if (dvar < 1.0E-3)
            goto COMPLETE;
        else
            goto ONE;
    }
    /*---------------------------------------------------------------------------
    * If RCI_request=3, then apply the preconditioner on the vector
    * tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1]
    *---------------------------------------------------------------------------
    * NOTE that ipar[21] and ipar[22] contain FORTRAN style addresses,
    * therefore, in C code it is required to subtract 1 from them to get C style
    * addresses
    *---------------------------------------------------------------------------*/
    if (RCI_request == 3)
    {
        if (ipar[3] == 3)
        {
            tmp[ipar[22] - 1 + 0] = -2.0;
            tmp[ipar[22] - 1 + 1] = 0.08519601586107672;
            tmp[ipar[22] - 1 + 2] = -1.1590871369607090;
            tmp[ipar[22] - 1 + 3] = -0.65791939687456980;
            tmp[ipar[22] - 1 + 4] = 0.75660051476696133;
        }
        else
        {
            if (ipar[3] == 4)
            {
                tmp[ipar[22] - 1 + 0] = 0.0;
                tmp[ipar[22] - 1 + 1] = 0.0;
                tmp[ipar[22] - 1 + 2] = 0.0;
                tmp[ipar[22] - 1 + 3] = 1.0;
                tmp[ipar[22] - 1 + 4] = -1.0;
            }
            else
            {
                for (i = 0; i < N; i++)
                {
                    tmp[ipar[22] - 1 + i] = (double) (i) * tmp[ipar[21] - 1 + i];
                }
            }
        }
        goto ONE;
    }
    /*---------------------------------------------------------------------------
    * If RCI_request=4, then check if the norm of the next generated vector is
    * not zero up to rounding and computational errors. The norm is contained
    * in dpar[6] parameter
    *---------------------------------------------------------------------------*/
    if (RCI_request == 4)
    {
        if (dpar[6] < 1.0E-12)
            goto COMPLETE;
        else
            goto ONE;
    }
    /*---------------------------------------------------------------------------
    * If RCI_request=anything else, then dfgmres subroutine failed
    * to compute the solution vector: computed_solution[N]
    *---------------------------------------------------------------------------*/
    else
    {
        goto FAILED;
    }
    /*---------------------------------------------------------------------------
    * Reverse Communication ends here
    * Get the current iteration number and the FGMRES solution (DO NOT FORGET to
    * call dfgmres_get routine as computed_solution is still containing
    * the initial guess!). Request to dfgmres_get to put the solution
    * into vector computed_solution[N] via ipar[12]
    *---------------------------------------------------------------------------*/
COMPLETE:
    ipar[12] = 0;
    dfgmres_get (&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
    /*---------------------------------------------------------------------------
    * Print solution vector: computed_solution[N] and the number of iterations: itercount
    *---------------------------------------------------------------------------*/
    printf (" The system has been solved \n");
    printf ("\n The following solution has been obtained: \n");
    for (i = 0; i < N; i++)
    {
        printf ("computed_solution[%d]=", i);
        printf ("%e\n", computed_solution[i]);
    }
    printf ("\n The expected solution is: \n");
    for (i = 0; i < N; i++)
    {
        printf ("expected_solution[%d]=", i);
        printf ("%e\n", expected_solution[i]);
        expected_solution[i] -= computed_solution[i];
    }
    printf ("\n Number of iterations: %d\n", itercount);
    i = 1;
    dvar = dnrm2 (&ivar, expected_solution, &i);

    /*-------------------------------------------------------------------------*/
    /* Release internal MKL memory that might be used for computations         */
    /* NOTE: It is important to call the routine below to avoid memory leaks   */
    /* unless you disable MKL Memory Manager                                   */
    /*-------------------------------------------------------------------------*/
    MKL_Free_Buffers ();

    if (itercount == expected_itercount && dvar < 1.0e-14)
    {
        printf ("\nThis example has successfully PASSED through all steps of ");
        printf ("computation!\n");
        return 0;
    }
    else
    {
        printf ("\nThis example may have FAILED as either the number of iterations differs\n");
        printf ("from the expected number of iterations %d, ", expected_itercount);
        printf ("or the computed solution\n");
        printf ("differs much from the expected solution (Euclidean norm is %e), or both.\n", dvar);
        return 1;
    }
    /*-------------------------------------------------------------------------*/
    /* Release internal MKL memory that might be used for computations         */
    /* NOTE: It is important to call the routine below to avoid memory leaks   */
    /* unless you disable MKL Memory Manager                                   */
    /*-------------------------------------------------------------------------*/
FAILED:
    printf ("\nThis example FAILED as the solver has returned the ERROR code %d", RCI_request);
    MKL_Free_Buffers ();
    #endif
    std::cout<<"Stub function, not implemented yet\n";
    return true;
}
int pardiso_solver::solve_sys_bcg(linear_sys &sys)
{
    std::cout<<"Stub function, not implemented yet\n";
    return false;
}
