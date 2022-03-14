#ifndef ECS_CLIENT_NOT_INTERSECT_H
#define ECS_CLIENT_NOT_INTERSECT_H


// Dependencies: find_std_vector(), vector_select()

auto not_intersect(Eigen::VectorXd m,Eigen::VectorXd n){
    // m - (N_mx1) column vector
    // n - (N_nx1) column vector

    // size of input vectors
    int N_m = m.size(); // number of entries in vector m
    int N_n = n.size(); // number of entries in vector n

    // variable declaration
    int N_l, N_s;
    Eigen::VectorXd* lvec_ptr; // pointer to larger vector
    Eigen::VectorXd* svec_ptr; // pointer to smaller vector   
    
    if (N_m > N_n){
        N_l = N_m;
        N_s = N_n;
        lvec_ptr = &m;
        svec_ptr = &n;
    }
    else{
        N_l = N_n;
        N_s = N_m;
        lvec_ptr = &n;
        svec_ptr = &m;
    }

    // variable declaration
    std::vector<int> index_l; // index of intersecting entries in larger vector
    std::vector<int> index_s; // index of intersecting entries in smaller vector

    // Find index of overlapping entries
    for (int j = 0; j < N_l; ++j){
        auto lvec_jth = (*lvec_ptr)(j); // j-th entry in larger vector
        for (int k = 0; k < N_s; ++k){
            auto svec_kth = (*svec_ptr)(k); // k-th entry in larger vector
            if (svec_kth == lvec_jth){
                index_l.push_back(j);
                index_s.push_back(k);
            }
        }
    }

    // Create vectors of Ones. (flag = 1).
    std::vector<int> flag_l(N_l,1); // flags for larger vector
    std::vector<int> flag_s(N_s,1); // flags smaller vector

    // set intersecting/overlapping entries to 0
    for (int i: index_l){
        flag_l.at(i) = 0;
    }
    for (int i: index_s){
        flag_s.at(i) = 0;
    }  

    // indices of non-overlapping entries
    std::vector<int> idx_nonov_l = find_std_vector(flag_l); // index of non-overlapping entries of larger vector
    std::vector<int> idx_nonov_s = find_std_vector(flag_s); // index of non-overlapping entries of smaller vector

    Eigen::VectorXd lvec = *lvec_ptr; // larger vector
    Eigen::VectorXd svec = *svec_ptr; // smaller vector

    // non-intersecting entries
    Eigen::VectorXd nonov_l = vector_select(lvec,idx_nonov_l); // non-overlapping entries of larger vector
    Eigen::VectorXd nonov_s = vector_select(svec,idx_nonov_s); // non-overlapping entries of smaller vector


    std::vector<int> index_nonov_m; // index of non-intersecting entries of vector m
    std::vector<int> index_nonov_n; // index of non-intersecting entries of vector n

    if (N_m > N_n){
        index_nonov_m = idx_nonov_l;
        index_nonov_n = idx_nonov_s;
    }
    else{
        index_nonov_m = idx_nonov_s;
        index_nonov_n = idx_nonov_l;               
    }

    // output struct declaration
    struct returnVals{
        std::vector<int> idx_nonov_m;
        std::vector<int> idx_nonov_n;
    };
    return returnVals{index_nonov_m,index_nonov_n};
}
#endif // ECS_CLIENT_NOT_INTERSECT_H