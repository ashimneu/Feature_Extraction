#ifndef ECS_CLIENT_SAMPLE_POINTS_H
#define ECS_CLIENT_SAMPLE_POINTS_H


std::vector<PointType> sample_points (){

    // create vector of points to contain pcd points
    std::vector<PointType> vector_of_points;
    // vector_of_points.reserve(10);

    Eigen::Vector4f point_L1(1,1,1,1);
    size_t laser_idx1 = 0;
    float range1 = 10.0;
    float hor_ang_d1 = 10.0;

    Eigen::Vector4f point_L2(2,2,2,1);
    size_t laser_idx2 = 1;
    float range2 = 20.0;
    float hor_ang_d2 = 20.0;

    Eigen::Vector4f point_L3(3,3,3,1);
    size_t laser_idx3 = 2;
    float range3 = 30.0;
    float hor_ang_d3 = 30.0;

    Eigen::Vector4f point_L4(4,4,4,1);
    size_t laser_idx4 = 0;
    float range4 = 11.0;
    float hor_ang_d4 = 11.0;

    Eigen::Vector4f point_L5(5,5,5,1);
    size_t laser_idx5 = 1;
    float range5 = 21.0;
    float hor_ang_d5 = 21.0;

    Eigen::Vector4f point_L6(6,6,6,1);
    size_t laser_idx6 = 2;
    float range6 = 31.0;
    float hor_ang_d6 = 31.0;

    Eigen::Vector4f point_L7(7,7,7,1);
    size_t laser_idx7 = 0;
    float range7 = 12.0;
    float hor_ang_d7 = 12.0;

    Eigen::Vector4f point_L8(8,8,8,1);
    size_t laser_idx8 = 1;
    float range8 = 22.0;
    float hor_ang_d8 = 22.0;

    Eigen::Vector4f point_L9(9,9,9,1);
    size_t laser_idx9 = 2;
    float range9 = 32.0;
    float hor_ang_d9 = 32.0;

    // push into vector of points
    vector_of_points.emplace_back(point_L1, laser_idx1, range1, hor_ang_d1);
    vector_of_points.emplace_back(point_L2, laser_idx2, range2, hor_ang_d2);
    vector_of_points.emplace_back(point_L3, laser_idx3, range3, hor_ang_d3);

    vector_of_points.emplace_back(point_L4, laser_idx4, range4, hor_ang_d4);
    vector_of_points.emplace_back(point_L5, laser_idx5, range5, hor_ang_d5);
    vector_of_points.emplace_back(point_L6, laser_idx6, range6, hor_ang_d6);

    vector_of_points.emplace_back(point_L7, laser_idx7, range7, hor_ang_d7);
    vector_of_points.emplace_back(point_L8, laser_idx8, range8, hor_ang_d8);
    vector_of_points.emplace_back(point_L9, laser_idx9, range9, hor_ang_d9);

    return vector_of_points;
}

#endif // ECS_CLIENT_SAMPLE_POINTS_H