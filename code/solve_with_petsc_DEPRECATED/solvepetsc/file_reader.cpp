#include<iostream>
#include "file_reader.h"


/*###########################
file_reader class implementation:
#############################*/

file_reader::file_reader(const std::string filename):
i_f(filename)
{
    if(!i_f.is_open())
    {
        std::cout<<"ERR in file reader: can't read file specified\n";
        //throw; //error handling not a priority at this time, using throw might supress potential optimizations
    }
}

file_reader::~file_reader()
{
    i_f.close();
}

bool file_reader::is_mat_asym()
{
    if(order_check++!=0)
    {
        std::cout<<"ERR: file reading functions not called in correct order(is_mat_asym)";
        //throw;  //might supress optimizations
    }
    bool nasym;
    i_f>>nasym;
    return nasym;
}

size_t file_reader::get_mat_sz()
{
    if(order_check++!=1)
    {
        std::cout<<"ERR: file reading functions not called in correct order(get_mat_sz)\n";
        std::cout<<"order_check is "<<order_check<<std::endl;
        //throw;  //might supress optimizations
    }
    size_t dummy;
    i_f>>dummy;
    return dummy;
}
size_t file_reader::get_non_diag_no()
{
    if(order_check++!=2)
    {
        std::cout<<"ERR: file reading functions not called in correct order(get_non_diag_no)";
        //throw;  //commented out as might supress optimizations
    }
    size_t non_diag_no;
    i_f>>non_diag_no;
    return non_diag_no;
}


std::deque<size_t> file_reader::read_row_index(size_t non_diag_no)
{
    size_t temp;
    std::deque<size_t> vect;

    if(order_check++!=3)
    {
        std::cout<<"ERR: file reading functions not called in correct order(1)";
        //throw;  //commented out as might supress optimizations
    }
    //as we know the size of the vector we can preallocate vecto dimension:
    vect.reserve(non_diag_no);
    for(size_t i = 0; i<non_diag_no && i_f.good();i++)
    {
        i_f>>temp;
        vect.push_back(temp);
    }
    if(!i_f.good())
    {
        std::cout<<"ERR in file reader, read row index\n";
    }
    return vect;
}
std::deque<size_t> file_reader::read_col_change(size_t sz )
{
    size_t temp;
    std::deque<size_t> vect;
    if(order_check++!=4)
    {
        std::cout<<"ERR: file reading functions not called in correct order(2)";
        //throw;  //commented out as might supress optimizations
    }
    //as we know the size of the vector we can preallocate vecto dimension:
    vect.reserve(sz+1);
    for(size_t i = 0; i<sz+1 && i_f.good();i++)
    {
        i_f>>temp;
        vect.push_back(temp);
    }
    if(!i_f.good())
    {
        std::cout<<"ERR in file reader, read col change\n";
    }
    return vect;
}
std::deque<double> file_reader::read_diag(size_t sz)
{
    double temp;
    std::deque<double> vect;

    if(order_check++!=5)
    {
        std::cout<<"ERR: file reading functions not called in correct order(3)";
        //throw;  //commented out as might supress optimizations
    }
    vect.reserve(sz);
    for(size_t i = 0; i<sz && i_f.good();i++)
    {
        i_f>>temp;
        vect.push_back(temp);
    }
    if(!i_f.good())
    {
        std::cout<<"ERR in file reader, read diag elem\n";
    }
    return vect;
}
std::deque<double> file_reader::read_non_diag(size_t non_diag_no,bool nasym)
{
    double temp;
    std::deque<double> vect;

    if(order_check++!=6)
    {
        std::cout<<"ERR: file reading functions not called in correct order(4)";
        //throw;  //commented out as might supress optimizations
    }
    if(!nasym)
    {
        vect.reserve(non_diag_no);
        for(size_t i = 0; i<non_diag_no && i_f.good();i++)
        {
            i_f>>temp;
            vect.push_back(temp);

        }
        if(!i_f.good())
        {
            std::cout<<"ERR in file reader, read non-diag elem\n";
        }
    }
    else
    {
        size_t temp2 = non_diag_no*2;
        vect.reserve(temp2);
        for(size_t i = 0; i<temp2 && i_f.good();i++)
        {
            i_f>>temp;
            vect.push_back(temp);
        }
        if(!i_f.good())
        {
            std::cout<<"ERR in file reader, read non-diag elem\n";
        }
    }
    return vect; //should be using move semantics so this is efficient
}
std::vector<double> file_reader::read_rhs(size_t sz)
{
    double temp;
    std::vector<double> vect;
    if(order_check++!=7)
    {
        std::cout<<"ERR: file reading functions not called in correct order(5)";
        //throw;  //commented out as might supress optimizations
    }
    vect.reserve(sz);
    for(size_t i = 0; i<sz && i_f.good();i++)
    {
        i_f>>temp;
        vect.push_back(temp);
    }
    //no checks of i_f good here as it is ambiguous wether the input file contains a solution or not. this will be checked in file_reader::read_sol
    return vect; //should be using move semantics so this is efficient
}
std::vector<double> file_reader::read_sol(size_t mat_dim)
{
    double temp;
    std::vector<double> vect;

    if(order_check++!=8)
    {
        std::cout<<"ERR: file reading functions not called in correct order(6)";
        //throw;  //commented out as might supress optimizations
    }
    if(i_f.good())
    {
        vect.reserve(mat_dim);
        for(size_t i = 0; i<mat_dim && i_f.good();i++)
        {
            i_f>>temp;
            vect.push_back(temp);
        }


        i_f>>temp; //see if there is still one more elem or this is the end
        if(i_f.good())
        {
            std::cout<<"ERR: file should not contain any more elements at this point..\n";
        }
    }
    else
    {
        std::cout<<"warning:seems as no solution array provided in input file(or err in read nsz vector)\n";
    }
    return vect;

}


/*#############################################*/
