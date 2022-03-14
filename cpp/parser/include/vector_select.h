#ifndef ECS_CLIENT_VECTOR_SELECT_H
#define ECS_CLIENT_VECTOR_SELECT_H

// extract entries from a vector of Eigen::VectorXd class.

auto vector_select(Eigen::VectorXd input_vector, std::vector <int> selection_index){
    
    int num_sel = selection_index.size(); // number of entries to be selected
    Eigen::VectorXd output_vector (num_sel);

    for (int i = 0; i < num_sel; ++i){
        output_vector(i) = input_vector(selection_index.at(i));    
    }
    return output_vector;
}
#endif // ECS_CLIENT_VECTOR_SELECT_H