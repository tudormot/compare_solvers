#ifndef LIN_SYS_H_INCLUDED
#define LIN_SYS_H_INCLUDED

class sparse_mat{
/*a class representing a matrix. This class can't be instantiated, only it's derivatives */
public:
    size_t mat_dim;
    bool is_asym;
    size_t nnz; //no of non zeros of this spare mat;
protected:
    sparse_mat(); //makes it so this can't be instantiated
};


class mat_CCSdiag:public sparse_mat{
/* class representing a matrix in Compressed Column Storage with extracted diagonal format (as provided in input)*/
public:

    //Decided to use double ended queues rather than vectors for the CCSdiag format. This formati is only used temporarily before a conversion to
    //CRS. The input matrices can be very big (10 GB), as a result having the matrix in both formats even temporarily (during the conversion process)
    //should be a problem even for nodes with a lot of RAM memory. Although deques will be created slower than vectors, the advatage is that
    //they can be shrinked efficiently as the vectors of the other format increase, hence the program should have a smaller memory footprint this way

    std::deque<size_t> row_i;
    std::deque<size_t> col_ch;
    std::deque<double> non_diag;
    std::deque<double> diag;

    mat_CCSdiag(file_reader fread); //constructs matrix in this format from file
    mat_CCSdiag() = delete;
};

class mat_CRS:public sparse_mat{
/* class representing a matrix in Compressed Row Storage (as required by PETSc)*/
private:
    void CCSdiag_to_CRS(std::deque<size_t> &row_i,
                    std::deque<size_t> &col_ch,
                    std::deque<double> &non_diag,
                    std::deque<double> &diag, //vectors storing input (all of above)
                    std::vector<int> &col_i,
                    std::vector<int> &row_ch,
                    std::vector<double> &elem_new, //vectors storing output
                    const bool is_asym); /*used in constructor*/

public:
    std::vector<size_t> row_i;
    std::vector<size_t> col_ch;
    std::vector<double> non_diag;
    std::vector<double> diag;

    mat_CRS(mat_CCSdiag&& old_mat); //constructor used to creat a mat in CRS format from a mat in CCSdiag format . Only r-values accepted for the time being
    mat_CRS() = delete;
};




#endif // LIN_SYS_H_INCLUDED
