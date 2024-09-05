#include <iostream>

#include <Eigen/Dense>
#include "math/math.h"

using namespace math;

bool floatEquals(double a, double b)
{
    return std::abs(a - b) < 1e-6;
}

double deg_to_rad = EIGEN_PI / 180.0;

double rad_to_deg = 57.2957795;


Eigen::Matrix3d rotate_x(double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           1.0, 0.0, 0.0,
            0.0, std::cos(radians), -std::sin(radians),
            0.0, std::sin(radians), std::cos(radians);
    return matrix;
}


Eigen::Matrix3d rotate_y(double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           std::cos(radians), 0.0, std::sin(radians),
            0.0, 1.0, 0.0,
            -std::sin(radians), 0.0, std::cos(radians);
    return matrix;
}


Eigen::Matrix3d rotate_z(double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           std::cos(radians), -std::sin(radians), 0.0,
            std::sin(radians), std::cos(radians), 0.0,
            0.0, 0.0, 1.0;
    return matrix;
}
Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    matrix = I * rotate_z(e(0)) * rotate_y(e(1)) * rotate_x(e(2));
    return matrix;
}
Eigen::Vector3d euler_zyx_from_rotation(const Eigen::Matrix3d &r)
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    if(floatEquals(r(2,0), -1.0))
    {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0,1), r(1,1));
    }
    else if(floatEquals(r(2,0), 1.0))
    {
        b = -(EIGEN_PI / 2.0);
        a = 0.0;
        c = -std::atan2(r(0,1), r(1,1));
    }
    else
    {
        b = std::atan2(-r(2,0), std::sqrt(r(0,0) * r(0,0) + r(1,0) * r(1,0)));
        a = std::atan2(r(1,0), r(0,0));
        c = std::atan2(r(2,1), r(2,2));
    }
    return Eigen::Vector3d(a, b, c);
}




int main()
{
    Eigen::Vector3d e = Eigen::Vector3d(60.0, 45.0, 30.0) * deg_to_rad;
    Eigen::Matrix3d rotation_matrix= math::rotation_matrix_from_euler_zyx(e);
    Eigen::Vector3d ea = euler_zyx_from_rotation(rotation_matrix);
    std::cout << " E:" << e.transpose() * rad_to_deg << std::endl;
    std::cout << "Ea:" << ea.transpose() * rad_to_deg << std::endl;
    return 0;
}
