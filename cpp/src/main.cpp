
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <string>
#include <typeinfo>
#include <Eigen/Dense>
#include <vector>

#include "arguparser.h"
#include "line_parser.h"
#include "point_type.h"
#include "clustering.h"

#define MAX_LINE_SIZE 500000
#define ESTIMATED_NUM_OF_POINTS 50000
#define MAX_BLOCK_SIZE 20
#define DELIMITER ','

int
main(int argc, char const * argv[])
{
    using namespace feature_detect;

    // retrieve pointcloud file pathing through command argument
    commandline_args cmd_args;
    if(program_opt_parser(argc, argv, cmd_args))
    {
        return EXIT_SUCCESS;
    }
    std::cout << "pointcloud file path: " << cmd_args.pc_file_path << std::endl;

    // load pcd
    pcl::PointCloud<pcl::PointXYZ>::Ptr input_pointcloud(new pcl::PointCloud<pcl::PointXYZ>);
    if(pcl::io::loadPCDFile<pcl::PointXYZ>(cmd_args.pc_file_path, *input_pointcloud) == -1)
    {
        std::cerr << "failed to read pointcloud file. \n";
        return -1;
    }
    else
        std::cout << "file loaded, size: " << (*input_pointcloud).size() << "\n";

    // extract scanline number from file_path
    std::string scanline_number_string = cmd_args.pc_file_path;
    scanline_number_string.erase(0, 3); // erase first three characters - "sl_"
    scanline_number_string.erase(scanline_number_string.size() - 4,
                                 4); // erase last four characters - ".pcd"
    int scanline_num{ std::stoi(scanline_number_string) };

    // load logID_log_P.txt
    std::ifstream log_file("laserID_log_sl.txt");
    char one_scanline_log[MAX_LINE_SIZE];
    char block_data[MAX_BLOCK_SIZE];
    size_t current_position{ 0 };

    while(log_file.getline(one_scanline_log, MAX_LINE_SIZE)) // get new line and check scanline num
    {
        current_position = 0; // reset and parse first block from line to get scanline number
        current_position = parse_block(one_scanline_log, current_position, block_data,
                                       MAX_BLOCK_SIZE, DELIMITER);
        if (std::stoi(block_data) >= scanline_num) // keep parsing until the sl_num is read
            break;
    }

    // create vector of points to contain all pcd points
    std::vector<PointType> vector_of_points;
    vector_of_points.reserve(ESTIMATED_NUM_OF_POINTS);

    if(std::stoi(block_data) == scanline_num)
    {
        size_t point_index{ 0 }; // start with the first point
        std::cout << "correct header found, scanline number = " << scanline_num << "\n";
        while(one_scanline_log[current_position]) // current position of scanline log is valid
        {
            // get point in L-frame from pointcloud
            auto & pc_point = input_pointcloud->points[point_index];
            Eigen::Vector4f point_L(pc_point.x, pc_point.y, pc_point.z, 1);

            // get corresponding log information
            current_position = parse_block(one_scanline_log, current_position, block_data,
                                           MAX_BLOCK_SIZE, DELIMITER);
            size_t laser_idx = std::stoi(block_data);
            current_position = parse_block(one_scanline_log, current_position, block_data,
                                           MAX_BLOCK_SIZE, DELIMITER);
            float range = std::stof(block_data);
            current_position = parse_block(one_scanline_log, current_position, block_data,
                                           MAX_BLOCK_SIZE, DELIMITER);
            float hor_ang_d = std::stof(block_data) * 57.29578f;

            // push into vector of points
            vector_of_points.emplace_back(point_L, laser_idx, range, hor_ang_d);
            ++point_index;
        }
        /*
        // test outputs
        std::cout << "scanline: " << scanline_num << "\n";
        for(auto const & point: vector_of_points)
        {
            std::cout << point.p_L.transpose() << "," << point.laser_idx << ", " << point.range
                      << ", " << point.hor_ang_d << "\n";
        }
         */

        /*    
        // convert std::vector to Eigen::Vector
        auto* ptr_vec_pnt = &vector_of_points;
        size_t num_pnt = vector_of_points.size();
        Eigen::Map<Eigen::VectorXd> vector_of_points_EC(ptr_vec_pnt, num_pnt); // acronym: EC - Eigen Class
        */
    }
    else
    {
        std::cout << "error, scanline number " << scanline_number_string << " not found in log. \n";
    }

    // Test codes
    Eigen::MatrixXd m1(2,3);
    m1 << 1,2,3,4,5,6;
    std::cout << m1(1,1) << std::endl

    // do processing using pcd - input: pointcloud, output: indices of corners

    // extract corner points using indices, push into new pointcloud of corners
    pcl::PointCloud<pcl::PointXYZ> corner_points; // container for new corner points
    // sample insert: corner_points.emplace_back(corner_point);

    // save corner points as a new pcd file
    if(!corner_points.empty())
    {
        pcl::io::savePCDFileASCII("corner_points.pcd", corner_points);
        std::cout << "corner points successfully saved. size: " << corner_points.size() << "\n";
    }
    return 0;
}
