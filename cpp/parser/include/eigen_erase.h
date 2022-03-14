#ifndef ECS_CLIENT_EIGEN_ERASE_H
#define ECS_CLIENT_EIGEN_ERASE_H

Eigen::MatrixXd eigen_erase(Eigen::MatrixXd input_vector, std::vector<int> idx_erase){
    // 1 dimensional matrix representing a vector (dims allowed: 1xN or Nx1)
    // idx_erase - indices of entries to be removed from vector.


    int rows  = input_vector.rows();
    int cols  = input_vector.cols();
    int num_v = input_vector.size(); // number of entries in input_vector
    int num_i = idx_erase.size();    // number of indices
    int dnum  = num_v - num_i;

    std::vector<int> idx_all = sequence(0,num_v-1,1);
    auto [v]

    if (rows > cols){
        Eigen::MatrixXd output_vector(dnum,1);


    }
    else if (rows < cols){
        Eigen::MatrixXd output_vector(1,dnum);

    } 
    else{

    }


    std::vector<int> sequence;

    if (increment > 0 ){
        for (int i = start; i <= end; i +=increment){
            sequence.push_back(i);
        }
    }
    else{
        for (int i = start; i >= end; i += increment){
            sequence.push_back(i);
        }
    }
    return sequence;
}
#endif // ECS_CLIENT_EIGEN_ERASE_H