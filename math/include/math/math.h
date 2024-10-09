#ifndef MATH_H
#define MATH_H

#include "Eigen/Dense"
#include "cstring"

namespace math{
    constexpr double deg_to_rad = EIGEN_PI / 180.0;
    constexpr double rad_to_deg = 180.0 / EIGEN_PI;
    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);
    Eigen::Vector3d skew_symmetric_to_vector(const Eigen::Matrix3d &m);
    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                    const Eigen::Vector3d &y,
                                                    const Eigen::Vector3d &z);
    Eigen::Matrix3d rotate_x(double radians);
    Eigen::Matrix3d rotate_y(double radians);
    Eigen::Matrix3d rotate_z(double radians);
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians);
    Eigen::Matrix3d rotation_matrix_from_euler(const char *s, const Eigen::Vector3d &radians);
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);
    void transform_vector();
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta);
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);
    bool floatEquals(double a, double b);
    Eigen::Vector3d euler_zyx_from_rotation(const Eigen::Matrix3d &r);
    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);
    double cot(double radians);
    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t);
}

#endif
