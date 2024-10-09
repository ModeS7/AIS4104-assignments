#include <iostream>

#include <Eigen/Dense>
#include <cstring>

#include "math/math.h"

// Most important rule: if its theta its in degrees.

Eigen::Matrix3d math::skew_symmetric(Eigen::Vector3d v)
{
    // Turn 3d vector into skew symmetric matrix. Modern Robotics page 75.
    Eigen::Matrix3d skew_matrix;
    skew_matrix <<
                0.0, -v(2), v(1),
            v(2), 0.0, -v(0),
            -v(1), v(0), 0.0;
    return skew_matrix;
}

Eigen::Vector3d math::skew_symmetric_to_vector(const Eigen::Matrix3d &m)
{
    // Turn skew symmetric matrix back to 3d vector. Modern Robotics page 75.
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
    // Create a 3x3 rotation matrix around x-axis. Modern Robotics page 72.
    double c = std::cos(radians);
    double s = std::sin(radians);
    Eigen::Matrix3d matrix;
    matrix <<
           1.0, 0.0, 0.0,
            0.0, c, -s,
            0.0, s, c;
    return matrix;
}

Eigen::Matrix3d math::rotate_y(double radians)
{
    // Create a 3x3 rotation matrix around y-axis. Modern Robotics page 72.
    double c = std::cos(radians);
    double s = std::sin(radians);

    Eigen::Matrix3d matrix;
    matrix <<
    c, 0.0, s,
    0.0, 1.0, 0.0,
    -s, 0.0, c;
    return matrix;
}

Eigen::Matrix3d math::rotate_z(double radians)
{
    // Create a 3x3 rotation matrix around z-axis. Modern Robotics page 72.
    double c = std::cos(radians);
    double s = std::sin(radians);

    Eigen::Matrix3d matrix;
    matrix <<
    c, -s, 0.0,
    s, c, 0.0,
    0.0, 0.0, 1.0;
    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians)
{
    // Create a 3x3 rotation matrix from an axis and an angle. Modern Robotics page 72.
    double c = std::cos(radians);
    double s = std::sin(radians);

    Eigen::Matrix3d matrix;
    matrix <<
    c + axis(0)*axis(0)*(1-c),
    axis(0)*axis(1)*(1-c)-axis(2)*s,
    axis(0)*axis(2)*(1-c)+axis(1)*s,
    axis(1)*axis(0)*(1-c)+axis(2)*s,
    c + axis(1)*axis(1)*(1-c),
    axis(1)*axis(2)*(1-c)-axis(0)*s,
    axis(2)*axis(0)*(1-c)-axis(1)*s,
    axis(2)*axis(1)*(1-c)+axis(0)*s,
    c + axis(2)*axis(2)*(1-c);
    return matrix;
}

Eigen::Matrix3d math::rotation_matrix_from_euler(const char *s, const Eigen::Vector3d &radians)
{
    // Create a 3x3 rotation matrix from euler angles
    // and given sequence. Modern Robotics page 577.
    Eigen::Matrix3d matrix;
    int order_code = -1;
    if (std::strcmp(s, "zyx") == 0)
        order_code = 0;
    else if (std::strcmp(s, "yzx") == 0)
        order_code = 1;
    else if (std::strcmp(s, "xzy") == 0)
        order_code = 2;
    else if (std::strcmp(s, "xyz") == 0)
        order_code = 3;
    else if (std::strcmp(s, "yxz") == 0)
        order_code = 4;
    else if (std::strcmp(s, "zxy") == 0)
        order_code = 5;

    switch (order_code) {
        case 0: // "zyx"
            matrix = math::rotate_z(radians(0)) *
                     math::rotate_y(radians(1)) *
                     math::rotate_x(radians(2));
            break;
        case 1: // "yzx"
            matrix = math::rotate_y(radians(0)) *
                     math::rotate_z(radians(1)) *
                     math::rotate_x(radians(2));
            break;
        case 2: // "xzy"
            matrix = math::rotate_x(radians(0)) *
                     math::rotate_z(radians(1)) *
                     math::rotate_y(radians(2));
            break;
        case 3: // "xyz"
            matrix = math::rotate_x(radians(0)) *
                     math::rotate_y(radians(1)) *
                     math::rotate_z(radians(2));
            break;
        case 4: // "yxz"
            matrix = math::rotate_y(radians(0)) *
                     math::rotate_x(radians(1)) *
                     math::rotate_z(radians(2));
            break;
        case 5: // "zxy"
            matrix = math::rotate_z(radians(0)) *
                     math::rotate_x(radians(1)) *
                     math::rotate_y(radians(2));
            break;
        default:
            // Invalid input: return identity matrix (or handle as necessary)
            matrix = Eigen::Matrix3d::Identity();
            break;
    }
    return matrix;
}

Eigen::Matrix4d math::transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
{
    // Create a 4x4 homogeneous transformation matrix
    // from a rotation matrix and a position vector.
    // Modern Robotics page 87.
    Eigen::Matrix4d matrix;
    matrix <<
    r(0, 0), r(0, 1), r(0, 2), p(0),
    r(1, 0), r(1, 1), r(1, 2), p(1),
    r(2, 0), r(2, 1), r(2, 2), p(2),
    0.0, 0.0, 0.0, 1.0;
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
// Rodrigues' formula for matrix exponential.
// Creates a rotation matrix from rotation axis
// w_hat and angle theta. Modern Robotics page 82.
{
    double rad = theta;
    Eigen::Matrix3d w_sk = math::skew_symmetric(w);
    return (Eigen::Matrix3d::Identity() + w_sk * std::sin(rad) + w_sk * w_sk * (1 - std::cos(rad)));
}

Eigen::Matrix4d math::matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
// Creat a homogeneous transformation from skew axis
// components(w and v) and theta angle. Modern Robotics page 103.
{
    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    const double w_norm = w.norm();
    const double v_norm = v.norm();
    const double rad = theta;
    if (w_norm == 1)
    {
        const Eigen::Matrix3d w_sk = math::skew_symmetric(w);
        const Eigen::Matrix3d ewt = math::matrix_exponential(w, theta);
        Eigen::Matrix3d G = Eigen::Matrix3d::Identity() * rad
                + (1 - std::cos(rad)) * w_sk + (rad - std::sin(rad)) * w_sk * w_sk;
        Eigen::Vector3d Gv = G * v;
        T <<
        ewt(0, 0), ewt(0, 1), ewt(0, 2), Gv(0),
        ewt(1, 0), ewt(1, 1), ewt(1, 2), Gv(1),
        ewt(2, 0), ewt(2, 1), ewt(2, 2), Gv(2),
        0.0, 0.0, 0.0, 1.0;
    }
    else if (v_norm == 1 and w_norm == 0)
    {
        T <<
        1.0, 0.0, 0.0, v(0),
        0.0, 1.0, 0.0, v(1),
        0.0, 0.0, 1.0, v(2),
        0.0, 0.0, 0.0, 1.0;

    }
    else
    {
        std::cout << "Invalid input" << std::endl;
        return Eigen::Matrix4d::Zero();
    }
    return T;
}

Eigen::MatrixXd math::adjoint_matrix(const Eigen::Matrix4d &tf)
// Create the adjoint representation 6x6 matrix from
// a transformation matrix. Modern Robotics page 98.
{
    Eigen::Matrix3d R = tf.block<3,3>(0,0);
    Eigen::Vector3d p = tf.block<3,1>(0,3);
    const Eigen::Matrix3d p_skew = math::skew_symmetric(p);

    Eigen::Matrix<double, 6, 6> adjoint_matrix;
    adjoint_matrix <<
    R, Eigen::Matrix3d::Zero(),
    p_skew * R, R;
    return adjoint_matrix;
}

bool math::floatEquals(double a, double b)
// Compare two floating point numbers with a tolerance.
{
    return std::abs(a - b) < 1e-6;
}

Eigen::Vector3d math::euler_zyx_from_rotation(const Eigen::Matrix3d &r)
// Extract euler angles from a rotation matrix.
// Modern Robotics page 579.
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
// Print the position and euler angles from transformation matrix.
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
// Cotangent function.
{
    return std::cos(radians) / std::sin(radians);
}

std::pair<Eigen::Vector3d, double> math::matrix_logarithm(const Eigen::Matrix3d &r)
// Compute the matrix logarithm of a rotation matrix.
// Modern Robotics page 85-86.
{
    Eigen::Matrix3d w_sk;
    double theta;
    double rad;
    Eigen::Vector3d w;
    if (r.isApprox(Eigen::Matrix3d::Identity(), 1e-6))
    {
        w = Eigen::Vector3d::Zero();
        std::cout << "w_sk is undefined" << std::endl;
        rad = 0.0;
    }
    else if (r.trace() == -1)
    {
        Eigen::Vector3d w_ = Eigen::Vector3d(r(0, 2), r(1, 2), r(2, 2)+1);
        w = w_ / std::sqrt(2 * (1 + r(2, 2)));
        rad = EIGEN_PI;
    }
    else
    {
        rad = (std::acos((r.trace() - 1) / 2));
        w_sk = r - r.transpose() / (2 * std::sin(rad));
        w = math::skew_symmetric_to_vector(w_sk);
    }
    theta = rad;
    return std::make_pair(w, theta);
}

std::pair<Eigen::VectorXd, double> math::matrix_logarithm(const Eigen::Matrix4d &t)
// Compute the matrix logarithm of a transformation matrix.
// Modern Robotics page 104.
{
    Eigen::Matrix3d R = t.block<3, 3>(0, 0);
    Eigen::Vector3d p = t.block<3, 1>(0, 3);
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double theta;
    if (R.isApprox(Eigen::Matrix3d::Identity(), 1e-6))
    {
        w = Eigen::Vector3d::Zero();
        v = p/p.norm();
        theta = p.norm();
    }
    else
    {
        std::pair<Eigen::Vector3d, double> m = math::matrix_logarithm(R);
        w = m.first;
        theta = m.second;
        Eigen::Matrix3d w_sk = math::skew_symmetric(w);
        double rad = theta;
        Eigen::Matrix3d G_1 = Eigen::Matrix3d::Identity()/rad - w_sk / 2 + (1 / rad - math::cot(rad / 2) / 2) * w_sk * w_sk;
        v = G_1 * p;
    }
    return std::make_pair(twist(w, v), theta);
}

































