#ifndef ECS_CLIENT_POINT_TYPE_H
#define ECS_CLIENT_POINT_TYPE_H

struct PointType
{
    Eigen::Vector4f p_L; // L-Frame coordinate of the point (X, Y, Z, 1)
    int laser_idx;       // source of the point
    float range;         // the Lidar laser range corresponding to this point
    float hor_ang_d;     // the Lidar laser horizontal angle corresponding to this point, degrees


    PointType(const Eigen::Vector4f & point_L, int laser_idx, float range, float hor_ang_d)
        : p_L(point_L), laser_idx(laser_idx), range(range), hor_ang_d(hor_ang_d)
    {
    };
};

#endif // ECS_CLIENT_POINT_TYPE_H
