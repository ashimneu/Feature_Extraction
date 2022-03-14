#ifndef ECS_CLIENT_INTERSECT_H
#define ECS_CLIENT_INTERSECT_H

auto intersect(Eigen::MatrixXd m,Eigen::MatrixXd n){

    // size of input vectors
    int N_m = m.size(); // number of entries in vector m
    int N_n = n.size(); // number of entries in vector n

    // variable declaration
    int N_l, N_s;
    Eigen::MatrixXd* large_vec_ptr;
    Eigen::MatrixXd* small_vec_ptr;    
    
    if (N_m > N_n){
        N_l = N_m;
        N_s = N_n;
        large_vec_ptr = &m;
        small_vec_ptr = &n;
    }
    else{
        N_l = N_n;
        N_s = N_m;
        large_vec_ptr = &n;
        small_vec_ptr = &m;
    }

    // variable declaration
    std::vector<int> val_intersect; // intersecting entries
    std::vector<int> index_l; // index of intersecting entries in larger vector
    std::vector<int> index_s; // index of intersecting entries in smaller vector

    // Find intersecting entries in large vec & small vec
    for (int j = 0; j < N_l; ++j){
        auto lvec_entry = (*large_vec_ptr)(0,j);
        for (int k = 0; k < N_s; ++k){
            auto svec_entry = (*small_vec_ptr)(0,k);
            if (svec_entry == lvec_entry){
                val_intersect.push_back(lvec_entry); // intersecting value
                index_l.push_back(j);
                index_s.push_back(k);
            }
        }
    }

    // variable declaration    
    std::vector<int> index_m; // index of intersecting entries in vector m
    std::vector<int> index_n; // index of intersecting entries in vector n

    if (N_m > N_n){
        index_m = index_l;
        index_n = index_s;
    }
    else{
        index_m = index_s;
        index_n = index_l;        
    }

    // return struct declaration
    struct returnVals{
        std::vector<int> val_intersect;
        std::vector<int> index_m;
        std::vector<int> index_n;
    };

    return returnVals{val_intersect,index_m,index_n};
}
#endif // ECS_CLIENT_INTERSECT_H