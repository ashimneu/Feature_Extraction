
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <typeinfo>

#include "sequence.h"
#include "intersect.h"
#include "vector_select.h"
#include "matrix_select.h"
#include "find_std_vector.h"
#include "not_intersect.h"
#include "point_type.h"
#include "sample_points.h"

int
main()
{
    using namespace std;
    Eigen::MatrixXd test_matrix(5,3);
    test_matrix << 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;
    cout << test_matrix << endl;
    cout << "rows: " << test_matrix.rows() << endl;
    cout << "cols: " << test_matrix.cols() << endl;
    cout << "col type: " << typeid(test_matrix.cols()).name()<<endl;
    //auto p = (test_matrix.cols()).cast<int>();
    //cout << "casted col type: " << (typeid(p).name()).cast<<endl
    //cout << "size: " << test_matrix.size() << endl;

    //Eigen::MatrixXd test_matrix2(num_pnt,1);
    //test_matrix2 << 1,2,3,4,5;
    cout << "innerSize" << test_matrix.innerSize() << endl;
    cout << "outerSize" << test_matrix.outerSize() << endl;
    cout << "innerStride" << test_matrix.innerStride() << endl;
    cout << "outerStride" << test_matrix.outerStride() << endl;

    // threshold check
    Eigen::ArrayXXi flag3 = (test_matrix.array() > 3).cast<int>();
    cout << flag3 << endl;

    // test intersect function
    Eigen::MatrixXd m(1,5), n(1,8);
    m << 1,2,3,4,5;
    n << 3,4,5,6,7,8,9,10;
    auto [val_intersect,idx_fl_m,idx_fl_n] = intersect(m,n); // flagged entries & their indices.

    //std::vector<int> ind(1,2);

    //Eigen::VectorXd m2(10);
    Eigen::MatrixXd m2(3,5);
    m2 << 0,1,2,3,4,
          5,6,7,8,9,
          10,11,12,13,14;
    auto a = m2.block(0,1,3,1);
    cout << a << endl;

    //
    cout << "idx_nonov_m";
    std::vector<int> idx_nonov_m{sequence(0,10-1,1)};
    for (int i: idx_nonov_m){
        cout << i << ' ';
    }
    cout << endl;


    cout << "------ Vector Select TEST ------" << endl;
    
    std::vector<int> selector;
    selector.push_back(2);
    selector.push_back(3);
    selector.at(0) = 4; 
    for (auto i:selector){
        cout << i << " ";
    }
    cout << endl;
    
    
    Eigen::VectorXd test_vec(5); // non-overlapping entries of larger vector
    test_vec << 100,101,102,103,104;
    cout << vector_select(test_vec,selector) << endl;


    cout << "------ Matrix Select TEST ------" << endl;
    Eigen::MatrixXd test_mtx(4,4); // non-overlapping entries of larger vector
    test_mtx << 1,2,3,4,
                5,6,7,8,
                9,10,11,12,
                13,14,15,16;
    std::vector<int> row_idx;
    std::vector<int> col_idx;
    std::vector<int> row_idx2;
    std::vector<int> col_idx2;

    row_idx.push_back(0); 
    row_idx.push_back(2);
    col_idx.push_back(1); 
    col_idx.push_back(3);
    row_idx2.push_back(-1);
    col_idx2.push_back(-1);
    
    //cout << test_mtx << endl;
    auto selected_mtx = matrix_select(test_mtx,row_idx,col_idx2);
    cout << selected_mtx << endl;                
    cout << endl;


    // comparison
    cout << "------ Comparison TEST ------" << endl;
    std::vector<int> z;
    z.push_back(99);

    cout << (z.at(0) == 99) << endl;

    int qq = 3;
    bool compare = qq <= test_matrix.cols();
    cout << "comparison: " << compare << endl;

    // find_std_vector function
    cout << "------ Find TEST ------" << endl;
    std::vector<int> test_integers;
    test_integers.push_back(9);
    test_integers.push_back(0);
    test_integers.push_back(5);
    test_integers.push_back(-0);
    test_integers.push_back(-1);

    auto selected_text_integers = find_std_vector(test_integers);
    for (int i: selected_text_integers){
        cout << i << " ";
    }
    cout << endl;

    // find non-overlapping entries of two vectors
    cout << "------ Non Intersecting Entries TEST ------" << endl;
    Eigen::VectorXd m3(5), n3(8);
    m3 << 1,2,3,4,5;
    n3 << 3,4,5,6,7,8,9,10;
    auto [idx_m,idx_n] = not_intersect(m3,n3); // indices
    
    cout << "idx_m: ";
    for (int i: idx_m){
        cout << i << " ";
    }
    cout << endl;

    cout << "idx_n: ";
    for (int i: idx_n){
        cout << i << " ";
    }
    cout << endl;

    // find non-overlapping entries of two vectors
    cout << "------ PointType ------" << endl;
    // create vector of points to contain pcd points
    std::vector<PointType> pc_struct = sample_points();
    std::vector<int> lasernum;
    
    for( std::vector< Foo >::const_iterator it = pc_struct.begin(); it != pc_struct.end(); ++it )
    {
        lasernum.push_back(it->laser_idx);
    }

    print_pc_struct(pc_struct);

    ///////////////////////////////////////////////////////////////////////
    /*
    // boolean array
    typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> myMatrixXd;
    myMatrixXd M(1,1);
    M << 9;
    cout << M << endl;

    // boolean array
    typedef Eigen::Array<bool,Eigen::Dynamic,1> ArrayXb;
    
    // bool test 1
    ArrayXb flag1(5);
    flag1 << true,true,true,true,true;
    cout << flag1 << endl;

    // bool test 2
    Eigen::MatrixXd test_matrix2(5,1);
    test_matrix2 << 3,3,3,3,3;
    ArrayXb flag2(5);
    // flag2 << test_matrix <= test_matrix2;
    //cout << flag2 << endl;
    */  

    return 0;
}
