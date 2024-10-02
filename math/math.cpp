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

Eigen::Vector3d math::skew_symetric_to_vector(const Eigen::Matrix3d &m)
{
    Eigen::Vector3d v;
    v << m(2, 1), m(0, 2), m(1, 0);
    return v;
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

Eigen::VectorXd math::screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
// Create a screw axis from a point, a direction and pitch. Modern Robotics page 101.
{
    Eigen::VectorXd screw_axis(6);
    screw_axis << s, - s.cross(q) + h * s;
    return screw_axis;
}

Eigen::Matrix3d math::matrix_exponential(const Eigen::Vector3d &w, double theta)
{
    Eigen::Matrix3d matrix;
    double rad = theta * math::deg_to_rad;
    Eigen::Matrix3d w_sk = math::skew_symmetric(w);
    matrix = Eigen::Matrix3d::Identity() + w_sk * std::sin(rad) + w_sk * w_sk * (1 - std::cos(rad));
    return matrix;
}

Eigen::Matrix4d math::matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
// Creat a matrix from a screw.
{
    Eigen::Matrix4d T;
    double rad = theta * math::deg_to_rad;
    if (w.norm() == 1)
    {
        Eigen::Matrix3d w_sk = math::skew_symmetric(w);
        Eigen::Matrix3d ewt = math::matrix_exponential(w, theta);
        Eigen::Matrix3d G = math::I_3() * rad + (1 - std::cos(rad)) * w_sk + (rad - std::sin(rad)) * w_sk * w_sk;
        T.block<3, 3>(0, 0) = ewt;
        T.block<3, 1>(0, 3) = G * v;
        T.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
        T(3, 3) = 1.0;
    }
    else if (v.norm() == 1 and w.norm() == 0)
    {
        T.block<3, 3>(0, 0) = math::I_3();
        T.block<3, 1>(0, 3) = v;
        T.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
        T(3, 3) = 1.0;
    }
    else
    {
        std::cout << "Invalid input" << std::endl;
        return Eigen::Matrix4d::Zero();
    }
    return T;
}