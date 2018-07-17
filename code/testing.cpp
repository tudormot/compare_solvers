#include<string>
#include<fstream>
#include<cmath>
#include<math.h>
#include<iostream>
#include<iomanip>
#include"sys_mat.h"
#include"testing.h"

#define DEBUG 0

bool test::test_input_file_dim(const linear_sys &sys, std::string filename)
{
    std::cout<<"Testing input file dimensions:\n";
    std::ifstream i_f;
    i_f.open(filename);
    i_f.clear();
    if(!i_f.is_open())
    {
     std::cout<<"ERR, file is not open\n";
    }
    size_t nasym, neq, nzs, expected_file_no_of_elem ,elem_read ;//number of equations, number of non zero elem
    double dummy;
    i_f>>nasym;
    i_f>>neq;
    i_f>>nzs;
    std::cout<<"File overview: nasym: "<<nasym<<"#neq: "<<neq<<"#nzs: "<<nzs<<'\n';
    if(nasym == 0)
    {
        expected_file_no_of_elem = 2* nzs + 3 *neq+1;

    }
    else
    {
        expected_file_no_of_elem = 3* nzs + 3 *neq+1;
    }
    std::cout<<"NOTE: this testing routing assumes that the input includes the solution vector as well\n";
    expected_file_no_of_elem=expected_file_no_of_elem+neq;
    elem_read = 0;
    while(i_f>>dummy && i_f.good())
    {
        elem_read ++;
    }

    std::cout<<"number of elements read from file: "<<elem_read<<'\n';
    std::cout<<"Expected no of elem to be read: "<<expected_file_no_of_elem<<'\n';
    if(elem_read == expected_file_no_of_elem)
    {
         std::cout<<"PASS\n";
         return true;
    }
    else
    {
        std::cout<<"FAIL\n";
        return false;
    }

}

bool test::analyze_col_ch(const linear_sys &sys)
/*analyze_col_ch: displays the last element of col_ch and checks if it is equal to neq+1, as the standard..*/
{
    std::cout<<"Last elem of col_ch is: "<<sys.col_ch.back()<<std::endl;
    //size_t temp = ;
    std::cout<<"We were expecting col_ch[end] to be "<<sys.non_diag_no+1<<std::endl;
    if(sys.non_diag_no +1== sys.col_ch.back())
    {
        std::cout<<"PASS\n";
        return true;
    }
    else
    {
        std::cout<<"FAIL\n";
        return false;
    }
}
bool test::check_row_i_range(const linear_sys &sys)
/* we would expect that in a symmetric matrix no row indexes below the diagonal are stored..this test checks that */
{
    bool symmetry  = true;
    bool is_upper = false;

    if(sys.col_ch[1]==sys.col_ch[2])//we need this check as we are using the second column to determine wether we expect the matrix to be upper or lower
    {
        std::cout<<"ERR in test::check_row_i_range: this particular matric has no elems on column 1: this test should be rewritten for this case as well\n";
        return false;
    }
    if(sys.row_i[sys.col_ch[1]]<2)
    {
        is_upper = true;
        std::cout<<"DEBUG: matrix is_upper = true\n";
    }
    else
    {
        std::cout<<"matrix has elems in lower triagle. If matrix is also symmetric, then it means elems are stored in lower triangle\n";
    }
    //not very efficient implementation but it's a test..
    std::cout<<"DEBUG: sys.mat_dim = "<<sys.mat_dim<<'\n';
    for(size_t i = 0;i<sys.mat_dim-1;i++)
    {
        for(size_t j = sys.col_ch[i]-1;j < sys.col_ch[i+1]-1;j++) //-1 stems rom the fact that col_ch contains indexes in fortran format
        {
            if((is_upper && sys.row_i[j]>i) || (!is_upper && sys.row_i[j]<i+1))
            {
                std::cout<<"DEBUG: col_ch[i],col_ch[i+1] = "<<sys.col_ch[i]<<' '<<sys.col_ch[i+1]<<'\n';
                std::cout<<"DEBUG: i,j,sys.row_i[j] = "<<i<<' '<<j<<' '<<sys.row_i[j]<<'\n';
                symmetry = false;
                break;
            }
        }
    }
    if(!symmetry != sys.is_asymmetric)
    {
        std::cout<<" test detected mismatch between matrix format and is_assymetric boolean input. sys.is_asymmetric = "<<sys.is_asymmetric<<std::endl;
        std::cout<<"FAIL\n";
        return false;
    }
    else
    {
        std::cout<<"PASS\n";
        return true;
    }


    return true;


}
bool test::calculate_sol_tolerance(double* true_sol, double* calc_sol, int vec_size)
   /* calculates the relative tolerance as the 2norm of solr-solc/2norm of solr*/
{
    double rel_tol = 0;
    double true_sol_norm= 0;
    //this could be written much better:
    for(int i =0;i<vec_size;i++)
    {
        rel_tol=rel_tol+fabs(true_sol[i]-calc_sol[i])/static_cast<double>(vec_size);
        true_sol_norm=true_sol_norm+fabs(true_sol[i])/static_cast<double>(vec_size);
    }
    if(std::isinf(true_sol_norm))
    {
        std::cout<<"problem in calculation of relative error in solution: norm of solution vector is inf\n";
    }
    rel_tol=rel_tol/true_sol_norm;
    std::cout<<"relative tolerance of solution is: "<<rel_tol<<"\n";

    //test returns true no matter what, the tester should read the console output for this test
    return true;
}
bool test::check_row_i_no_diag(const linear_sys &sys)
{
    std::cout<<"warning, test not implemented, returning false\n";
    return false;
}

bool test::create_structsymmat_from_symmat(linear_sys &sys)
{
	std::ofstream myfile("my_structsym_mat.txt");
	if(!myfile.is_open())
	{
		std::cout<<"ERR in print to file. Can't create/open file specified\n";
		throw; //error handling not a priority at this time
	}
	if(sys.is_asymmetric)
	{
		std::cout<<"ERR: you are using an already structural symmetric matrix as input.symmetric matrix required\n";
		throw;
	}
	myfile<<!sys.is_asymmetric<<' '<<sys.mat_dim<<' '<<sys.non_diag_no<<' ';

	for(int i = 0;i<sys.non_diag_no;i++)
	{
		myfile <<sys.row_i[i]<<' ';
	}
	for(int i = 0;i<sys.mat_dim+1;i++)
	{
		myfile <<sys.col_ch[i]<<' ';
	}
	for(int i = 0;i<sys.mat_dim;i++)
	{
		myfile << std::fixed << std::setprecision(15) <<sys.diag[i]<<' ';
	}
	//add to file the lower triangular section:
	for(int i = 0;i<sys.non_diag_no;i++)
	{
		myfile << std::fixed << std::setprecision(15)  <<sys.non_diag[i]<<' ';
	}
	/*now to not make the matrix perfectly symmetric change some elems in non_diag vector*/
	//add some magic numbers to the matrix:
	//we know that input is bigger than 35000 as our smallest matrix is bigger than that
	for(int i = 1000; i <35000; i = i+1000)
	{
		sys.non_diag[i]=i/100;
	}
	//add to file upper triangular section:
	for(int i = 0;i<sys.non_diag_no;i++)
	{
		myfile<< std::fixed << std::setprecision(15)  <<sys.non_diag[i]<<' ';
	}
	//add rhs to file
	for(int i = 0;i<sys.mat_dim;i++)
	{
		myfile<< std::fixed << std::setprecision(15)  <<sys.rhs[i]<<' ';
	}
	//add sol to file
	for(int i = 0;i<sys.mat_dim;i++)
	{
		myfile<< std::fixed << std::setprecision(15)  <<sys.sol[i]<<' ';
	}
	myfile<<"\n";

	myfile.close();
	return true;
}


bool test::compare_matrices(linear_sys &sys1, linear_sys &sys2)
{
	if(sys1.mat_dim != sys2.mat_dim || sys1.non_diag_no != sys2.non_diag_no)
	{
		std::cout<<"ERROR: matrix dimensions do not agree\n";
		return false;
	}

	for(int i = 0;i<sys1.mat_dim;i++)
	{
		/*compare sols. (assumed that the input contained solution)*/
		if(sys1.sol[i]!=sys2.sol[i])
		{
			std::cout<<"ERROR:matrix solutions do not agree at index: "<<i<<std::endl;
			std::cout<<"sys1.sol is: "<<sys1.sol[i]<<" while sys2.sol is: "<<sys2.sol[i]<<std::endl;
			//return false;
		}
		/*compare rhs */
		if(sys1.rhs[i]!=sys2.rhs[i])
		{
			std::cout<<"ERROR:matrix rhs do not agree at index: "<<i<<std::endl;
			return false;
		}
		/*compare diag */
		if(sys1.diag[i]!=sys2.diag[i])
		{
			std::cout<<"ERROR:matrix diags do not agree at index: "<<i<<std::endl;
			return false;
		}
		/*compare col_ch*/
		if(sys1.col_ch[i]!=sys2.col_ch[i])
		{
			std::cout<<"ERROR:matrix col_ch do not agree at index: "<<i<<std::endl;
			return false;
		}
	}
	/*compare last entry of col_ch as well:*/
	if(sys1.col_ch[sys1.mat_dim] != sys2.col_ch[sys2.mat_dim])
	{
		std::cout<<"ERROR:matrix col_ch do not agree at index: "<<sys1.mat_dim<<std::endl;
		return false;
	}

	/*now compare row_i and elems in lower triangular matrix*/
	for(int i = 0;i<sys1.non_diag_no;i++)
	{
		if(sys1.row_i[i]!=sys2.row_i[i])
		{
			std::cout<<"ERROR: row i do not agree at index: "<<i<<std::endl;
			return false;
		}

		if(sys1.non_diag[i]!=sys2.non_diag[i])
		{
			std::cout<<"ERROR: non diag do not agree at index (both in lower triangular matrix): "<<i<<std::endl;
			return false;
		}
	}
	if(sys1.is_asymmetric == false && sys2.is_asymmetric == true)
	{
		for(int i = 0;i<sys1.non_diag_no;i++)
		{
			if(sys1.non_diag[i]!=sys2.non_diag[i+sys1.non_diag_no])
			{
				std::cout<<"ERROR: non diag do not agree at index (1st mat in lower triangular and 2nd mat in upper): "<<i<<std::endl;
				return false;
			}
		}
	}
	else
	{
		std::cout<<"ERROR! not implemented (in compare matrices test)\n";
	}
	return true;

}

void test::print_array_to_file(double * array,int sz,std::string output_file)
{
	//note: this should only be called from sequencial code, no check for node_ranks
	std::ofstream myfile(output_file);
    if(!myfile.is_open())
    {
        std::cout<<"ERR in print to file. Can't create/open file specified\n";
        throw; //error handling not a priority at this time
    }
	myfile<<sz<<' ';
	for(int i = 0;i<sz;i++)
	{
		myfile<< std::fixed << std::setprecision(15)  <<array[i]<<' ';
	}

	myfile.close();
}
