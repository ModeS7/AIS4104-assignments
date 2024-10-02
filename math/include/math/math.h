#ifndef MATH_H
#define MATH_H

#include "Eigen/Dense"
#include "cstring"

namespace math{
    constexpr double deg_to_rad = EIGEN_PI / 180.0;
    constexpr double rad_to_deg = 180.0 / EIGEN_PI;
    inline Eigen::Matrix3d I_3() {
        return Eigen::Matrix3d::Identity();
    }
    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);
    Eigen::Vector3d skew_symetric_to_vector(const Eigen::Matrix3d &m);
    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                    const Eigen::Vector3d &y,
                                                    const Eigen::Vector3d &z);
    Eigen::Matrix3d rotate_x(double radians);
    Eigen::Matrix3d rotate_y(double radians);
    Eigen::Matrix3d rotate_z(double radians);
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians);
    Eigen::Matrix3d rotation_matrix_from_euler(const char *s, const Eigen::Vector3d &e);
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p);
    void transform_vector();
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta);
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);
}

#endif
