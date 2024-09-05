#include <iostream>

#include <Eigen/Dense>
#include "math/math.h"

using namespace math;



bool floatEquals(double a, double b)
{
    return std::abs(a - b) < 1e-6;
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
    Eigen::Vector3d e = Eigen::Vector3d(60.0, 45.0, 30.0) * math::deg_to_rad;
    Eigen::Matrix3d rotation_matrix= math::rotation_matrix_from_euler_zyx(e);
    Eigen::Vector3d ea = euler_zyx_from_rotation(rotation_matrix);
    std::cout << " E:" << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Ea:" << ea.transpose() * math::rad_to_deg << std::endl;
    return 0;
}
