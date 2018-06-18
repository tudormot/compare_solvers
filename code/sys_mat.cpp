#include<iostream>
#include "sys_mat.h"

#define DEBUG (0)



linear_sys::linear_sys(const std::string filename,PetscInt total_no_nodes, PetscInt local_node_rank):
    no_of_nodes(total_no_nodes),
    node_rank(local_node_rank)
{
    if(node_rank==0)
    {
        /*warning: order matters, do not change*/
        file_reader f_in(filename);
        is_asymmetric=f_in.is_mat_asym();
        mat_dim=f_in.get_mat_sz();
        non_diag_no=f_in.get_non_diag_no();
        row_i=f_in.read_row_index(non_diag_no);
        col_ch=f_in.read_col_change(mat_dim);
        diag=f_in.read_diag(mat_dim);
        non_diag=f_in.read_non_diag(non_diag_no,is_asymmetric);
        rhs=f_in.read_rhs(mat_dim);
        sol=f_in.read_sol(mat_dim);

        //now broadcast the matrix sizes to the other linear_sys object (which won't contain the vectors, just the sizes)
        int packet_send[3];
        packet_send[0]=static_cast<int>(is_asymmetric);
        packet_send[1]=mat_dim;
        packet_send[2]=non_diag_no;
        MPI_Bcast(
            packet_send,
            3,
            MPI_INT,
            0,
            MPI_COMM_WORLD);

    }
    else
    {
        int packet_recv[2];
        MPI_Bcast(
            packet_recv,
            3,
            MPI_INT,
            0,
            MPI_COMM_WORLD);

        is_asymmetric=static_cast<bool>(packet_recv[0]);
        mat_dim=packet_recv[1];
        non_diag_no=packet_recv[2];

    }
}

linear_sys::~linear_sys()
{
#if DEBUG
    std::cout<<"dummy destruct\n";
#endif // DEBUG
}

void linear_sys::display_system_stats() const
{
    std::cout<<"Size of matrix (square mat) :"<<mat_dim<<'\n';
    std::cout<<"Symmetry of matrix (NOT sure what this means yet, just disaplying nasym) :"<<is_asymmetric<<'\n';
    std::cout<<"Number of non-zero non-diagonal elem :"<<non_diag_no<<'\n';
}

void linear_sys::release_mem_mat()
{
    if(this->node_rank==0)
    {
        row_i.resize(0);
        col_ch.resize(0);
        diag.resize(0);
        non_diag.resize(0);
        row_i.shrink_to_fit();
        col_ch.shrink_to_fit();
		diag.shrink_to_fit();
		non_diag.shrink_to_fit();

    }
}
//#################################################
//file reader implementation below
//#################################################
std::vector<int> file_reader::read_row_index(size_t non_diag_no)
{
    size_t temp;
    std::vector<int> vect;

    if(order_check++!=3)
    {
        std::cout<<"ERR: file reading functions not called in correct order(1)";
        throw;
    }
    //as we know the size of the vector we can preallocate vecto dimension:
    vect.reserve(non_diag_no);
    for(size_t i = 0; i<non_diag_no && i_f.good(); i++)
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
std::vector<int> file_reader::read_col_change(size_t sz )
{
    size_t temp;
    std::vector<int> vect;
    if(order_check++!=4)
    {
        std::cout<<"ERR: file reading functions not called in correct order(2)";
        throw;
    }
    //as we know the size of the vector we can preallocate vector dimension:
    vect.reserve(sz+1);
    for(size_t i = 0; i<sz+1 && i_f.good(); i++)
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
std::vector<double> file_reader::read_diag(size_t sz)
{
    double temp;
    std::vector<double> vect;

    if(order_check++!=5)
    {
        std::cout<<"ERR: file reading functions not called in correct order(3)";
        throw;
    }
    vect.reserve(sz);
    for(size_t i = 0; i<sz && i_f.good(); i++)
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
std::vector<double> file_reader::read_non_diag(size_t non_diag_no,bool nasym)
{
    double temp;
    std::vector<double> vect;

    if(order_check++!=6)
    {
        std::cout<<"ERR: file reading functions not called in correct order(4)";
        throw;
    }
    if(!nasym)
    {
        vect.reserve(non_diag_no);
        for(size_t i = 0; i<non_diag_no && i_f.good(); i++)
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
        for(size_t i = 0; i<temp2 && i_f.good(); i++)
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
        throw;
    }
    vect.reserve(sz);
    for(size_t i = 0; i<sz && i_f.good(); i++)
    {
        i_f>>temp;
        vect.push_back(temp);
    }
    //no checks of i_f good here as it is ambiguous whether the input file contains a solution or not. this will be checked in file_reader::read_sol
    return vect; //should be using move semantics so this is efficient
}
std::vector<double> file_reader::read_sol(size_t mat_dim)
{
    double temp;
    std::vector<double> vect;

    if(order_check++!=8)
    {
        std::cout<<"ERR: file reading functions not called in correct order(6)";
        throw;
    }
    i_f>>temp;
    if(i_f.good())
    {
        vect.reserve(mat_dim);
        for(size_t i = 0; i<mat_dim && i_f.good(); i++)
        {

            vect.push_back(temp);
            i_f>>temp;
        }

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

file_reader::file_reader(const std::string filename):
i_f(filename)
{
    if(!i_f.is_open())
    {
        std::cout<<"ERR in file reader: can't read file specified\n";
        throw; //error handling not a priority at this time
    }
}

file_reader::~file_reader()
{
    i_f.close();
}
int file_reader::get_mat_sz()
{
    if(order_check++!=1)
    {
        std::cout<<"ERR: file reading functions not called in correct order(get_mat_sz)\n";
        std::cout<<"order_check is "<<order_check<<std::endl;
        throw;
    }
    int dummy;
    i_f>>dummy;
    return dummy;
}
int file_reader::get_non_diag_no()
{
    if(order_check++!=2)
    {
        std::cout<<"ERR: file reading functions not called in correct order(get_non_diag_no)";
        throw;
    }
    int non_diag_no;
    i_f>>non_diag_no;
    return non_diag_no;
}
bool file_reader::is_mat_asym()
{
    if(order_check++!=0)
    {
        std::cout<<"ERR: file reading functions not called in correct order(is_mat_sym)";
        throw;
    }
    bool nasym;
    i_f>>nasym;
    return nasym;
}


