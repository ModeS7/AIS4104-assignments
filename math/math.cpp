#include <iostream>

#include <Eigen/Dense>
#include <cstring>

#include "math/math.h"

// Most important rule: if its theta its in degrees.

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
// Create a screw axis from a point, a direction and pitch.
// Where q is any point on the axis, s is a unit vector in the direction of the axis,
// and h is the pitch of the axis. Modern Robotics page 101.
{
    Eigen::VectorXd screw_axis(6);
    screw_axis << s, - s.cross(q) + h * s;
    return screw_axis;
}

Eigen::VectorXd math::twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
// Create a twist from two vectors: angular velocity w and linear velocity v.
// Modern Robotics page 96.
{
    Eigen::VectorXd twist(6);
    twist << w, v;
    return twist;
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

Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf)
// Create the adjoint matrix from a transformation matrix.
{
    Eigen::Matrix3d R = tf.block<3,3>(0,0);
    Eigen::Vector3d p = tf.block<3,1>(0,3);
    Eigen::MatrixXd adjoint_matrix(6,6);
    adjoint_matrix << R, Eigen::Matrix3d::Zero(),
            math::skew_symmetric(p) * R, R;
    return adjoint_matrix;
}

bool math::floatEquals(double a, double b)
{
    return std::abs(a - b) < 1e-6;
}

Eigen::Vector3d math::euler_zyx_from_rotation(const Eigen::Matrix3d &r)
{
    double a = 0.0;
    double b = 0.0;
    double c = 0.0;

    if(math::floatEquals(r(2,0), -1.0))
    {
        b = EIGEN_PI / 2.0;
        a = 0.0;
        c = std::atan2(r(0,1), r(1,1));
    }
    else if(math::floatEquals(r(2,0), 1.0))
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

void math::print_pose(const std::string &label, const Eigen::Matrix4d &tf)
{
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);
    Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
    Eigen::Vector3d e = math::euler_zyx_from_rotation(r);
    std::cout << label << std::endl;
    //std::cout << "Euler ZYX(rad): " << e.transpose() << std::endl;
    std::cout << "Euler ZYX(deg): " << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Position: " << p.transpose() << std::endl;
}

double math::cot(double radians)
{
    return std::cos(radians) / std::sin(radians);
}

std::pair<Eigen::Vector3d, double> math::matrix_logarithm(const Eigen::Matrix3d &r)
{
    Eigen::Matrix3d w_sk;
    double theta;
    Eigen::Vector3d w;
    if (r.isApprox(Eigen::Matrix3d::Identity(), 1e-6))
    {
        w_sk = Eigen::Matrix3d::Zero();
        std::cout << "w_sk is undefined" << std::endl;
        theta = 0.0;
    }
    else if (r.trace() == -1)
    {
        Eigen::Vector3d w_ = Eigen::Vector3d(r(0, 2), r(1, 2), r(2, 2)+1);
        w /= std::sqrt(2 * (1 + r(2, 2)));
        theta = 180.0;
    }
    else
    {
        theta = std::acos((r.trace() - 1) / 2);
        w_sk = r - r.transpose() / (2 * std::sin(theta * math::rad_to_deg));
        w = math::skew_symetric_to_vector(w_sk);
    }
    return std::make_pair(w, theta);
}

std::pair<Eigen::VectorXd, double> math::matrix_logarithm(const Eigen::Matrix4d &t)
{
    Eigen::Matrix3d R = t.block<3, 3>(0, 0);
    Eigen::Vector3d p = t.block<3, 1>(0, 3);
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double theta;
    if (R.isApprox(math::I_3(), 1e-6))
    {
        w = Eigen::Vector3d::Zero();
        v = p/p.norm();
        theta = p.norm() * rad_to_deg;
    }
    else
    {
        std::pair<Eigen::Vector3d, double> m = math::matrix_logarithm(R);
        w = m.first;
        theta = m.second;
        Eigen::Matrix3d w_sk = math::skew_symmetric(w);
        double rad = theta * math::deg_to_rad;
        Eigen::Matrix3d G_1 = math::I_3()/rad - w_sk / 2 + (1 / rad - math::cot(rad / 2) / 2) * w_sk * w_sk;
        v = G_1 * p;
        w = math::skew_symetric_to_vector(w_sk);
    }
    return std::make_pair(twist(w, v), theta);
}

































