#ifndef PARDISO_SOLVER_H_INCLUDED
#define PARDISO_SOLVER_H_INCLUDED

/*class that encapsulates the intel pardiso library*/
class pardiso_solver:public lin_sys_solver
{
private:

	/*vectors to store linear system in format accepted by Intel Pardiso:*/

	std::vector<int>col_i;  //this vector tells us where the col is being changed(in terms of indexes of ja)
    std::vector<int>row_ch;  //these are the row indexes of element.
    std::vector<double>elem;   //this vector contains all the non-zero values in the matrix


public:
    /*  virtual int solve_sys(linear_sys &sys):
     * solves system with the intel pardiso direct solver, also displays quite a lot of informative output
     * at the end of the solving process, checks whether the input file contained a solution.
     * If a pre-existing solution is found stored inside the linear_sys object, solve_sys calls
     * test::calculate_sol_tolerance to calculate and display the relative error between the solutions*/
    virtual int solve_sys(linear_sys &sys);

    int solve_sys_gmres(linear_sys &sys);   //stub function, as I have realised intel pardiso does not provide parallel gmres implementation
    int solve_sys_bcg(linear_sys &sys);     //stub function, as I have realised intel pardiso does not provide parallel bcg implementation
    virtual ~pardiso_solver();
    pardiso_solver(const linear_sys &sys); //constructor that also performs the conversion from original input to input format required by pardiso library
    virtual void print_sol_to_file(linear_sys &sys, std::string & outputfile);
};



#endif // PARDISO_SOLVER_H_INCLUDED
