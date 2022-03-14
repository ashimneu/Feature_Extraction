#ifndef ECS_CLIENT_FIND_H
#define ECS_CLIENT_FIND_H


// Find index of non-zero entries in a vector of std::vector<int> class.
// Dependencies: sequence()

auto find_std_vector(std::vector<int> input_vector){
    
    std::vector<int> output_vector; // declare output
    
    int num = input_vector.size(); // number of entries
    
    // sequence of integers from 0 to length(input_vector)
    std::vector<int> seq = sequence(0,num-1,1);

    for (int i = 0; i < num; ++i){

        //check if ith entry is non-zero
        if (input_vector.at(i) != 0){
            output_vector.push_back(seq.at(i));
        };
    }
    return output_vector;
}
#endif // ECS_CLIENT_FIND_H