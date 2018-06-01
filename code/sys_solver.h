#ifndef SYS_SOLVER_H_INCLUDED
#define SYS_SOLVER_H_INCLUDED

class lin_sys_solver{
protected:
     void CLCXformats_to_CRS(const std::vector<int> &row_i_old,const std::vector<int> &col_ch_old,const std::vector<double> &non_diag_old, const std::vector<double> &diag_old, //vectors storing input (variable name contains _old)
                              std::vector<int> &col_i_new,           std::vector<int> &row_ch_new     ,      std::vector<double> &elem_new, const bool is_asym);
     void print_array_to_file(double *array, int sz, std::string output_file);
     std::string output_file;
public:

    std::vector<double>sol;

    virtual int solve_sys(linear_sys &sys) = 0;
    virtual void print_sol_to_file(linear_sys &sys) = 0;
    virtual ~lin_sys_solver();


};

#include"intel_pardiso_solver.h"
#include"petsc_solver.h"


#endif // SYS_SOLVER_H_INCLUDED
