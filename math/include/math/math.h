//
// Created by ASAP on 05/09/2024.
//

#ifndef MATH_H
#define MATH_H

#include "Eigen/Dense"

namespace math{


    double deg_to_rad(double degrees);
    double rad_to_deg(double radians);
    Eigen::Matrix3d skew_symmetric(Eigen::Vector3d v);
    void skew_symmetric_test();
Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                    const Eigen::Vector3d &y,
                                                    const Eigen::Vector3d &z);
    Eigen::Matrix3d rotate_x(double degrees);

// c)
    Eigen::Matrix3d rotate_y(double degrees)
    {
        Eigen::Matrix3d matrix;
        double radians = deg_to_rad(degrees);
        matrix <<
               std::cos(radians), 0.0, std::sin(radians),
                0.0, 1.0, 0.0,
                -std::sin(radians), 0.0, std::cos(radians);
        return matrix;
    }

// d)
    Eigen::Matrix3d rotate_z(double degrees)
    {
        Eigen::Matrix3d matrix;
        double radians = deg_to_rad(degrees);
        matrix <<
               std::cos(radians), -std::sin(radians), 0.0,
                std::sin(radians), std::cos(radians), 0.0,
                0.0, 0.0, 1.0;
        return matrix;
    }

// e)
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees)
    {
        Eigen::Matrix3d matrix;
        double radians = deg_to_rad(degrees);
        matrix <<
               std::cos(radians)+axis(0)*axis(0)*(1-std::cos(radians)), axis(0)*axis(1)*(1-std::cos(radians))-axis(2)*std::sin(radians), axis(0)*axis(2)*(1-std::cos(radians))+axis(1)*std::sin(radians),
                axis(1)*axis(0)*(1-std::cos(radians))+axis(2)*std::sin(radians), std::cos(radians)+axis(1)*axis(1)*(1-std::cos(radians)), axis(1)*axis(2)*(1-std::cos(radians))-axis(0)*std::sin(radians),
                axis(2)*axis(0)*(1-std::cos(radians))-axis(1)*std::sin(radians), axis(2)*axis(1)*(1-std::cos(radians))+axis(0)*std::sin(radians), std::cos(radians)+axis(2)*axis(2)*(1-std::cos(radians));
        return matrix;
    }

// f)
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
    {
        Eigen::Matrix3d matrix;
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        matrix = I * rotate_z(e(0)) * rotate_y(e(1)) * rotate_x(e(2));
        return matrix;
    }

// g)
    void rotation_matrix_test()
    {
        Eigen::Matrix3d rot =
                rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
        Eigen::Matrix3d rot_aa =
                rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0, 0.5773503}, 120.0);
        Eigen::Matrix3d rot_fa =
                rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
                                                Eigen::Vector3d{-0.5, -0.5, 0.707107},
                                                Eigen::Vector3d{0.707107, -0.707107, 0.0});
        std::cout << "Rotation matrix from Euler: " << std::endl;
        std::cout << rot << std::endl << std::endl;
        std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
        std::cout << rot_aa << std::endl << std::endl;
        std::cout << "Rotation matrix from frame axes: " << std::endl;
        std::cout << rot_fa << std::endl << std::endl;
    }

// -----------------------------------------
// 2.3 Task 3
// a)
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
        Eigen::Matrix4d matrix;
        matrix.block<3, 3>(0, 0) = r;
        matrix.block<3, 1>(0, 3) = p;
        matrix.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
        matrix(3, 3) = 1.0;
        return matrix;
    }

// b)
    void transformation_matrix_test()
    {
        Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
        Eigen::Vector3d v{1.0, -2.0, 3.0};
        std::cout << "transformation_matrix: " << std::endl;
        std::cout << transformation_matrix(r, v) << std::endl;
    }

// c)
    void transform_vector()
    {
        Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{60, 45.0, 0.0});
        Eigen::Vector3d v{0.0, 0.0, 10.0};
        Eigen::Vector4d v_a{2.5, 3.0, -10.0, 1.0};
        Eigen::Vector4d W_cord = transformation_matrix(r, v) * v_a;
        std::cout << "W_coordinates: " << std::endl;
        std::cout << W_cord << std::endl;
    }

}

#endif
