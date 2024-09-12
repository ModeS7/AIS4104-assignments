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
    double rad = theta * math::deg_to_rad;
    Eigen::Matrix3d w_sk = math::skew_symmetric(w);
    matrix = Eigen::Matrix3d::Identity() + w_sk * std::sin(rad) + w_sk * w_sk * (1 - std::cos(rad));
    return matrix;
}

// b)
std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r)
{
    Eigen::Matrix3d w_sk;
    double theta;
    Eigen::Vector3d w;
    if (r.isApprox(Eigen::Matrix3d::Identity(), 1e-6))
    {
        w_sk = Eigen::Matrix3d::Zero(); // w_sk should be returned as undefined ask Aleksander
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

// c)
Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
{
    Eigen::Matrix4d T;
    double rad = theta * math::deg_to_rad;
    if (w.norm() == 1)
    {
        Eigen::Matrix3d w_sk = math::skew_symmetric(w);
        Eigen::Matrix3d ewt = matrix_exponential(w, theta);
        Eigen::Matrix3d G = math::I_3() * rad + (1 - std::cos(rad)) * w_sk + (rad - std::sin(rad)) * w_sk * w_sk;
        T.block<3, 3>(0, 0) = ewt;
        T.block<3, 1>(0, 3) = G * v;
        T.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
        T(3, 3) = 1.0;
    }
    else if (v.norm() == 1 and w.norm() == 0) // w = 0 not sure what that means ask Aleksander
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

// d)
std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t)
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
        theta = p.norm(); //not sure if i get this in deg or rad
    }
    else
    {
        std::pair<Eigen::Vector3d, double> m = matrix_logarithm(R);
        w = m.first;
        theta = m.second;
        Eigen::Matrix3d w_sk = math::skew_symmetric(w);
        double rad = theta * math::deg_to_rad;
        Eigen::Matrix3d G_1 = math::I_3()/rad - w_sk / 2 + (1 / rad - cot(rad / 2) / 2) * w_sk * w_sk;
        v = G_1 * p;
        w = math::skew_symetric_to_vector(w_sk);
    }
    return std::make_pair(twist(w, v), theta);
}

//-----------------------------------------------------------------------------------------------
// Task 4
// a)
void print_pose(const std::string &label, const Eigen::Matrix4d &tf)
{
    Eigen::Vector3d p = tf.block<3, 1>(0, 3);
    Eigen::Matrix3d r = tf.block<3, 3>(0, 0);
    Eigen::Vector3d e = euler_zyx_from_rotation(r);
    std::cout << label << std::endl;
    std::cout << "Euler ZYX: " << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Position: " << p.transpose() << std::endl;
}

// b)
Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions)
{
    double L1 = 10.0, L2 = 10.0, L3 = 10.0;
    Eigen::Matrix4d t01;
    t01.block<3, 3>(0, 0) = math::rotate_z(joint_positions[0] * math::deg_to_rad);
    t01.block<3, 1>(0, 3) = Eigen::Vector3d(0.0, 0.0, 0.0);
    t01.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t01(3, 3) = 1.0;
    Eigen::Matrix4d t12;
    t12.block<3, 3>(0, 0) = math::rotate_z(joint_positions[1] * math::deg_to_rad);
    t12.block<3, 1>(0, 3) = Eigen::Vector3d(L1, 0.0, 0.0);
    t12.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t12(3, 3) = 1.0;
    Eigen::Matrix4d t23;
    t23.block<3, 3>(0, 0) = math::rotate_z(joint_positions[2] * math::deg_to_rad);
    t23.block<3, 1>(0, 3) = Eigen::Vector3d(L2, 0.0, 0.0);
    t23.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t23(3, 3) = 1.0;
    Eigen::Matrix4d t34 = Eigen::Matrix4d::Identity();
    t34.block<3, 1>(0, 3) = Eigen::Vector3d(L3, 0.0, 0.0);
    return t01 * t12 * t23 * t34;
}

// c)
Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions)
{
    double L1 = 10.0, L2 = 10.0, L3 = 10.0;

}


int main()
{
    /*Eigen::Vector3d e = Eigen::Vector3d(60.0, 45.0, 30.0) * math::deg_to_rad;
    Eigen::Matrix3d rotation_matrix= math::rotation_matrix_from_euler("zyx", e);
    Eigen::Vector3d ea = euler_zyx_from_rotation(rotation_matrix);
    std::cout << " E:" << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Ea:" << ea.transpose() * math::rad_to_deg << std::endl;*/
    //task2_a();
    //task2_b();
    //print_pose("shit", Eigen::Matrix4d::Identity());
    //task4b
    /*std::vector<double> joint_positions = {10.0, -15.0, 2.75};
    Eigen::Matrix4d transform = planar_3r_fk_transform(joint_positions);
    print_pose("end position:", transform);*/
    return 0;
}
















