#ifndef FILE_READER_H_INCLUDED
#define FILE_READER_H_INCLUDED

/*helper object that encapsulates the file reading process..*/
class file_reader{
private:

    std::ifstream i_f;
    unsigned int order_check = 0; // as the chunks of the file need to be read in a particular order, the member functions of this object are expected to be called in a articular order. This variable helps check that this is indeed the case

public:

    size_t get_mat_sz();
    bool is_mat_asym();
    size_t get_non_diag_no();
    std::deque<size_t> read_row_index(size_t non_diag_nr);
    std::deque<size_t> read_col_change(size_t col_nr);
    std::deque<double> read_diag(size_t mat_dim);
    std::deque<double> read_non_diag(size_t non_diag_nr,bool is_asym);

    std::vector<double> read_rhs(size_t mat_dim);
    std::vector<double> read_sol(size_t mat_dim);


    file_reader(const std::string filename);
     ~file_reader();

public:
    file_reader() = delete;
    file_reader(const file_reader& dummy) = delete;
    const file_reader& operator=(const file_reader& dummy) = delete;
};

#endif // FILE_READER_H_INCLUDED
