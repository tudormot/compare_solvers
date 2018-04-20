#ifndef TESTING_H_INCLUDED
#define TESTING_H_INCLUDED

/*test consistencyof matrix, input file etc.*/
class test{
public:
   static bool test_input_file_dim(const linear_sys &sys, std::string filename);
   /* test_input_file_dim: checks if what was read from file actually is what we expect */
   static bool analyze_col_ch(const linear_sys &sys);
   /*analyze_col_ch: displays the last element of col_ch and checks if it is equal to neq+1, as the standard..*/
   static bool check_row_i_range(const linear_sys &sys);
   /* we would expect that in a symmetric matrix no row indexes below the diagonal are stored..this test checks that */
   static bool chech_row_i_no_diag(const linear_sys &sys);
   /* just a sanity check that row_i does not contain indexes which would point on items of the main diagonal */
};

#endif // TESTING_H_INCLUDED
