#include <iostream>

#include <Eigen/Dense>
#include "math/math.h"


using namespace math;


// Task 1
// a)
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

// b)
Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::VectorXd twist(6);
    twist << w, v;
    return twist;
}

// c)
Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    Eigen::VectorXd screw_axis(6);
    screw_axis << s, - s.cross(q) + h * s;
    return screw_axis;
}

// d)
Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d R = tf.block<3,3>(0,0);
    Eigen::Vector3d p = tf.block<3,1>(0,3);
    Eigen::MatrixXd adjoint_matrix(6,6);
    adjoint_matrix << R, Eigen::Matrix3d::Zero(),
                      math::skew_symmetric(p) * R, R;
    return adjoint_matrix;
}

// e)
double cot(double radians)
{
    return std::cos(radians) / std::sin(radians);
}

//-----------------------------------------------------------------------------------------------
// Task 2
// a)
void task2_a()
{
    Eigen::Vector3d f_w = Eigen::Vector3d(-30.0, 0.0, 0.0);
    Eigen::Vector3d m_s = Eigen::Vector3d(0.0, 0.0, 2.0);
    Eigen::Vector3d e_ws = Eigen::Vector3d(60.0, -60.0, 0.0) * math::deg_to_rad;
    Eigen::Vector3d e_sw = Eigen::Vector3d(-60.0, 60.0, 0.0);

    Eigen::Matrix3d R_ws = math::rotation_matrix_from_euler("yzx", e_ws);
    Eigen::Matrix3d R_sw = R_ws.transpose();

    Eigen::Vector3d f_s = R_sw * f_w;

    Eigen::Vector3d m_w = R_ws * m_s;

    std::cout << "f_w: " << f_w.transpose() << std::endl;
    std::cout << "t_w: " << m_w.transpose() << std::endl;
    std::cout << "f_s: " << f_s.transpose() << std::endl;
    std::cout << "m_s: " << m_s.transpose() << std::endl;
}

// b)
void task2_b()
{
    Eigen::Vector3d m = Eigen::Vector3d(0.0, 0.0, -0.75);
    Eigen::Vector3d f = Eigen::Vector3d(0.0, -6.0, 0.0);

    Eigen::VectorXd F_f = twist(m, f);

    std::cout << "F_f: " << F_f.transpose() << std::endl;
}

//-----------------------------------------------------------------------------------------------
// Task 3
// a)
Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta)
{
    Eigen::Matrix3d matrix;
    double wx = math::skew_symmetric(w)(0,0);
    double wy = math::skew_symmetric(w)(1,1);
    double wz = math::skew_symmetric(w)(2,2);
    double w_norm = w.norm();
    matrix = Eigen::Matrix3d::Identity() + (std::sin(w_norm * theta) / w_norm) * math::skew_symmetric(w) + ((1 - std::cos(w_norm * theta)) / (w_norm * w_norm)) * math::skew_symmetric(w) * math::skew_symmetric(w);
    return matrix;
}



int main()
{
    /*Eigen::Vector3d e = Eigen::Vector3d(60.0, 45.0, 30.0) * math::deg_to_rad;
    Eigen::Matrix3d rotation_matrix= math::rotation_matrix_from_euler("zyx", e);
    Eigen::Vector3d ea = euler_zyx_from_rotation(rotation_matrix);
    std::cout << " E:" << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Ea:" << ea.transpose() * math::rad_to_deg << std::endl;*/
    //task2_a();
    task2_b();
    return 0;
}
