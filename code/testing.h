#ifndef TESTING_H_INCLUDED
#define TESTING_H_INCLUDED

/*test consistencyof matrix, input file etc.*/
/*note, these should all be called from only one processor.*/
class test{
public:
   static bool test_input_file_dim(const linear_sys &sys, std::string filename);
   /* test_input_file_dim: checks if what was read from file actually is what we expect */
   static bool analyze_col_ch(const linear_sys &sys);
   /*analyze_col_ch: displays the last element of col_ch and checks if it is equal to neq+1, as the standard..*/
   static bool check_row_i_range(const linear_sys &sys);
   /* we would expect that in a symmetric matrix no row indexes below the diagonal are stored..this test checks that */
   static bool check_row_i_no_diag(const linear_sys &sys);
   /* just a sanity check that row_i does not contain indexes which would point on items of the main diagonal */
   static bool calculate_sol_tolerance(double* true_sol, double* calc_sol, int vec_size);
   /* calculates the relative tolerance as the 2norm of solr-solc/2norm of solr*/
   static bool create_structsymmat_from_symmat(linear_sys &sys);
   /*programatically creates a structurally symmetric matrix from the symmetric matrix stored in sys
    * . This should only be used once to create a small structural symmetric matrix on which I can test the solver
    * implementations*/
   static void print_array_to_file(double *array, int sz, std::string output_file);
   static bool compare_matrices(linear_sys &sys1, linear_sys &sys2);
};

#endif // TESTING_H_INCLUDED
