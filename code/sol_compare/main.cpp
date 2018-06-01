#include<iostream>
#include<fstream>
#include<cmath>
#include<math.h>

/* this program compares solutions provided in two input files.*/
/* format of input file is size of vector followed by all elems of vector*/

int main(int argc, char *argv[])
{

    if(argc != 3)
    {
        std::cout<<"invalid input";
    }

    std::ifstream i_f1(argv[1]);
    i_f1.clear();
    if(!i_f1.is_open())
    {
     std::cout<<"ERR, file 1 is not open\n";
        return 0;
    }
    std::ifstream i_f2(argv[2]);
    i_f2.clear();
    if(!i_f2.is_open())
    {
     std::cout<<"ERR, file 2 is not open\n";
        return 0;
    }
    size_t sz1,sz2;
    double scalar1,scalar2;
    i_f1>>sz1;
    i_f2>>sz2;
    if(sz1 != sz2)
    {
        std::cout<<"ERR, size of input files differs \n";  
        return 0;     
    }
    
    double rel_tol = 0;
    double true_sol_norm= 0;
    for(int i =0;i<sz1;i++)
    {   
        i_f1>>scalar1;
        i_f2>>scalar2;
        rel_tol=rel_tol+fabs(scalar1-scalar2)/static_cast<double>(sz1);
        true_sol_norm=true_sol_norm+fabs(scalar1)/static_cast<double>(sz1);
    }
    rel_tol=rel_tol/true_sol_norm;
    std::cout<<"relative tolerance of solution is: "<<rel_tol<<"\n";

    i_f1.close();
    i_f2.close();    

    return 0;
}
