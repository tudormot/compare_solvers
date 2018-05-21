#ifndef SYS_SOLVER_H_INCLUDED
#define SYS_SOLVER_H_INCLUDED

class lin_sys_solver{
protected:
     void CCSdiag_to_CRS(const std::vector<int> &row_i_old,const std::vector<int> &col_ch_old,const std::vector<double> &non_diag_old, const std::vector<double> &diag_old, //vectors storing input (variable name contains _old)
                              std::vector<int> &col_i_new,           std::vector<int> &row_ch_new     ,      std::vector<double> &elem_new, const bool is_asym);

public:

    std::vector<double>sol;

    virtual int solve_sys(linear_sys &sys) = 0;
    virtual ~lin_sys_solver();


};

#include"intel_pardiso_solver.h"
#include"petsc_solver.h"


#endif // SYS_SOLVER_H_INCLUDED
