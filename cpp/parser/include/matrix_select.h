#ifndef ECS_CLIENT_MATRIX_SELECT_H
#define ECS_CLIENT_MATRIX_SELECT_H

// Entries of matrix of Eigen::MatrixXd class are extracted
// using indices of selected rows & columns.
// if sel_rows_idx == -1, all rows are extracted.
// if sel_cols_idx == -1, all columns are extracted.

auto matrix_select(Eigen::MatrixXd input_matrix, std::vector <int> sel_rows_idx,std::vector <int> sel_cols_idx){
    
    int count_mtxrows = input_matrix.rows();
    int count_mtxcols = input_matrix.cols();
    int count_selrows = sel_rows_idx.size(); // number of selected rows
    int count_selcols = sel_cols_idx.size(); // number of selected cols

    if (sel_rows_idx.at(0) == -1){
        // extract all rows & selected column
        Eigen::MatrixXd output_matrix (count_mtxrows,count_selcols);        
        for (int i = 0; i < count_mtxrows; ++i){
            for (int j = 0; j < count_selcols; ++j){
                    output_matrix(i,j) = input_matrix(i,sel_cols_idx.at(j));
                }    
        }
        return output_matrix;
    }
    else if (sel_cols_idx.at(0) == -1){
        // extract selected rows & all columns
        Eigen::MatrixXd output_matrix (count_selrows,count_mtxcols);        
        for (int i = 0; i < count_selrows; ++i){
            for (int j = 0; j < count_mtxcols; ++j){
                    output_matrix(i,j) = input_matrix(sel_rows_idx.at(i),j);
                }    
        }
        return output_matrix;
    }
    else {
        // extract only selected rows & columns
        Eigen::MatrixXd output_matrix (count_selrows,count_selcols);        
        for (int i = 0; i < count_selrows; ++i){
            for (int j = 0; j < count_selcols; ++j){
                    output_matrix(i,j) = input_matrix(sel_rows_idx.at(i),sel_cols_idx.at(j));
                }    
        }
        return output_matrix;
    }
}
#endif // ECS_CLIENT_MATRIX_SELECT_H