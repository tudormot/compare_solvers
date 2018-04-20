#include<string>
#include<fstream>
#include<iostream>
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
    //TODO this is not efficient implementation..
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
bool chech_row_i_no_diag(const linear_sys &sys);
/* just a sanity check that row_i does not contain indexes which would point on items of the main diagonal */
