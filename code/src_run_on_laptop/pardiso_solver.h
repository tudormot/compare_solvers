#ifndef PARDISO_SOLVER_H_INCLUDED
#define PARDISO_SOLVER_H_INCLUDED

class pardiso_solver:public lin_sys_solver{
private:
    //note:pardiso requires these vectors to be of type int
    //this is the data (IE matrix and RHS), however stored in simple CRS(compressed sparse col) rather than CCS with extracted diagonal.
    std::vector<int>col_i;  //this vector tells us where the col is being changed(in termms of indexes of ja) SZ(ia)=mat_dim+1
    std::vector<int>row_ch;  //these are the row indexes of element. SZ(ja) = SZ(a)
    std::vector<double>elem;   //this vector contains all the non-zero values in the matrix


public:
    virtual int solve_sys(linear_sys &sys);
    virtual ~pardiso_solver();
    pardiso_solver(const linear_sys &sys);
};



#endif // PARDISO_SOLVER_H_INCLUDED
