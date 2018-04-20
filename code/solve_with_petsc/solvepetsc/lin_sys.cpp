
#include "file_reader.h"
#include "lin_sys.h"



mat_CCSdiag::mat_CCSdiag(file_reader fread):
is_asym(fread.is_mat_asym()),
mat_dim(fread.get_mat_sz()),
nnz(fread.get_non_diag_no()+mat_dim),
row_i(fread.read_row_index(nnz-mat_dim)), //nnz-mat_dim = number of non zero elem which are not on main diagonal
col_ch(fread.read_col_change(mat_dim)),
diag(fread.read_diag(mat_dim)),
non_diag(fread.read_non_diag(nnz-mat_dim,is_asym))
{}
