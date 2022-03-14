#ifndef ECS_CLIENT_CLUSTERING_H
#define ECS_CLIENT_CLUSTERING_H


# define THRESHOLD_RANGE // for range based cluster separation
# define THRESHOLD_ANGLE // for angle based cluster separation

// vector_of_points : Eigen::Vector4f point_L, int laser_idx, float range, float hor_ang_d

std::vector<int> clustering_with_range_theta(const std::vector<PointType>* ptr_vector_of_points)
{
    size_t num_pnt = *ptr_vector_of_points.size() -1; // total points in pointCloud

    
    // convert std::vector to Eigen::Vector
    Eigen::Map<Eigen::MatrixXd> points_xyz_idx_range_angle(*ptr_vector_of_points, num_pnt); // acronym: EC - Eigen Class
    


    Eigen::VectorXi<int> cluster_number(num_pnt);    // uninitialized variable to store cluster number for each point

    Eigen::MatrixXd idx = points_xyz_idx_range_angle.block(0,0,num_pnt,2);
    Eigen::MatrixXd points_xyz = points_xyz_idx_range_angle.block(0,0,num_pnt,2);    
    Eigen::MatrixXd range = points_xyz_idx_range_angle.block(0,0,num_pnt,2);
    Eigen::MatrixXd angle = points_xyz_idx_range_angle.block(0,0,num_pnt,2);

    // compute differences in range and angle values
    Eigen::MatrixXd diff_range = range.block(1,0,num_pnt,0) - range.block(0,0,num_pnt-1,0); // difference in range
    Eigen::MatrixXd diff_angle = angle.block(1,0,num_pnt,0) - angle.block(0,0,num_pnt-1,0); // difference in angle

    
    // initialize flag variables to store result of following threshold comparisons
    typedef Array<bool,Dynamic,1> ArrayXb;
    ArrayXb flag_range(num_pnts);
    ArrayXb flag_angle(num_pnts);

    // do comparison with threshold on range & angle differences
    ;

    for (int i = 0; i<= num_pnt-1; ++i)
    {
        double diffvec = 
    }
    std::vector<int> diff

 
 
    std::vector<int> cluster_number(num_pnt)
    return cluster_number
}

#endif // ECS_CLIENT_CLUSTERING_H
