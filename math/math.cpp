#include <iostream>

#include <Eigen/Dense>
#include <cstring>

#include "math/math.h"


Eigen::Matrix3d math::skew_symmetric(Eigen::Vector3d v)
{
    Eigen::Matrix3d skew_matrix;
    skew_matrix <<
                0.0, -v(2), v(1),
            v(2), 0.0, -v(0),
            -v(1), v(0), 0.0;
    return skew_matrix;
}

Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                const Eigen::Vector3d &y,
                                                const Eigen::Vector3d &z)
{
    Eigen::Matrix3d matrix;
    matrix << x, y, z;
    return matrix;
}

Eigen::Matrix3d math::rotate_x(double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           1.0, 0.0, 0.0,
            0.0, std::cos(radians), -std::sin(radians),
            0.0, std::sin(radians), std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_y(double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           std::cos(radians), 0.0, std::sin(radians),
            0.0, 1.0, 0.0,
            -std::sin(radians), 0.0, std::cos(radians);
    return matrix;
}

Eigen::Matrix3d math::rotate_z(double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           std::cos(radians), -std::sin(radians), 0.0,
            std::sin(radians), std::cos(radians), 0.0,
            0.0, 0.0, 1.0;
    return matrix;
}


Eigen::Matrix3d math::rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians)
{
    Eigen::Matrix3d matrix;
    matrix <<
           std::cos(radians)+axis(0)*axis(0)*(1-std::cos(radians)), axis(0)*axis(1)*(1-std::cos(radians))-axis(2)*std::sin(radians), axis(0)*axis(2)*(1-std::cos(radians))+axis(1)*std::sin(radians),
            axis(1)*axis(0)*(1-std::cos(radians))+axis(2)*std::sin(radians), std::cos(radians)+axis(1)*axis(1)*(1-std::cos(radians)), axis(1)*axis(2)*(1-std::cos(radians))-axis(0)*std::sin(radians),
            axis(2)*axis(0)*(1-std::cos(radians))-axis(1)*std::sin(radians), axis(2)*axis(1)*(1-std::cos(radians))+axis(0)*std::sin(radians), std::cos(radians)+axis(2)*axis(2)*(1-std::cos(radians));
    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_euler(const char *s, const Eigen::Vector3d &e)
{
    Eigen::Matrix3d matrix;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    if(std::strcmp(s, "zyx") == 0) {
        matrix = I * math::rotate_z(e(0)) * math::rotate_y(e(1)) * math::rotate_x(e(2));
    }
    else if(std::strcmp(s, "yzx") == 0) {
        matrix = I * math::rotate_y(e(0)) * math::rotate_z(e(1)) * math::rotate_x(e(2));
    }
    else if(std::strcmp(s, "xzy") == 0) {
        matrix = I * math::rotate_x(e(0)) * math::rotate_z(e(1)) * math::rotate_y(e(2));
    }
    else if(std::strcmp(s, "xyz") == 0) {
        matrix = I * math::rotate_x(e(0)) * math::rotate_y(e(1)) * math::rotate_z(e(2));
    }
    else if(std::strcmp(s, "yxz") == 0) {
        matrix = I * math::rotate_y(e(0)) * math::rotate_x(e(1)) * math::rotate_z(e(2));
    }
    else if(std::strcmp(s, "zxy") == 0) {
        matrix = I * math::rotate_z(e(0)) * math::rotate_x(e(1)) * math::rotate_y(e(2));
    }
    return matrix;
}

Eigen::Matrix4d math::transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    Eigen::Matrix4d matrix;
    matrix.block<3, 3>(0, 0) = r;
    matrix.block<3, 1>(0, 3) = p;
    matrix.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    matrix(3, 3) = 1.0;
    return matrix;
}

