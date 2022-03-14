#ifndef ECS_CLIENT_PRINT_PC_STRUCT_H
#define ECS_CLIENT_PRINT_PC_STRUCT_H

void print_pc_struct(std::vector<PointType> vector_of_points){
    
    for(auto const & point: vector_of_points)
    {
        std::cout << point.p_L.transpose() << "," << point.laser_idx << ", " << point.range
                  << ", " << point.hor_ang_d << std::endl;
    }
}
#endif // ECS_CLIENT_PRINT_PC_STRUCT_H