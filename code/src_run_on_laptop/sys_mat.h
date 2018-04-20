#ifndef SYS_MAT_H_INCLUDED
#define SYS_MAT_H_INCLUDED

#include<vector>
#include<string>
#include<fstream>
#include<petscksp.h> //required for MPI routines and typedefs



/*helper class trying to tidy up the construction of the sparse_mat class..object will be deleted after sparse_mat construction process so it doesnt have to be memory efficient :) */
//non parallel object which should only be called from 1 process: IE process with rank 0
class file_reader{
private:

    std::ifstream i_f;
    unsigned int order_check = 0; // as the chunks of the file need to be read in a particular order, the member functions of this object are expected to be called in a articular order. This variable helps check that this is indeed the case

public:
    int get_mat_sz();
    bool is_mat_asym();
    int get_non_diag_no();
    std::vector<int> read_row_index(size_t non_diag_nr);
    std::vector<int> read_col_change(size_t col_nr);
    std::vector<double> read_diag(size_t mat_dim);
    std::vector<double> read_non_diag(size_t non_diag_nr,bool is_asym);
    std::vector<double> read_rhs(size_t mat_dim);
    std::vector<double> read_sol(size_t mat_dim);


    file_reader(const std::string filename);
    virtual ~file_reader();

public:
    file_reader() = delete;
    file_reader(const file_reader& dummy) = delete;
    const file_reader& operator=(const file_reader& dummy) = delete;
};

//a parallel object which contains different things depending on which process it is being created
//if it is being created on process rank 0, it will contain the input vectors and the sparse matrix sizes
//if it is being created on any other process, it will only contain the sparse matrix sizes(and some empty vectors TODO not very efficient) ( required for setting up PETSc)
class linear_sys{
private:

public:

    //TODO should these be public?
    bool is_asymmetric = false; //true if matrix is asymmetric
    int mat_dim; //matrices always square
    int non_diag_no;
    std::vector<int>row_i;//row index
    std::vector<int>col_ch;//column change
    std::vector<double>diag; //diag elements not part of the CCR
    std::vector<double>non_diag;
    std::vector<double>rhs;
    std::vector<double> sol; //this might be already filled with solution or not, depending on input textfile

    PetscInt    no_of_nodes, node_rank;

    linear_sys(const std::string filename, PetscInt total_no_nodes, PetscInt local_node_rank);
    linear_sys() = delete;
    linear_sys(const linear_sys& dummy) = delete;
    const linear_sys& operator=(const linear_sys& dummy) = delete;
    virtual ~linear_sys();

    //todo "parallelise" these functions
    void release_mem_mat(); //releases memory from the vectors that are stoing the matrix in "pardiso format"
    void display_system_stats() const;

};





#endif // SYS_MAT_H_INCLUDED
