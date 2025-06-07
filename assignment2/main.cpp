#include <iostream>

#include <Eigen/Dense>
#include "math/math.h"
#include <fstream>


using namespace math;


// Task 1
// a)

/*
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
        w = math::skew_symmetric_to_vector(w_sk);
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
        theta = p.norm() * rad_to_deg;
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
    std::cout << "Euler ZYX(rad): " << e.transpose() << std::endl;
    std::cout << "Euler ZYX(deg): " << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Position: " << p.transpose() << std::endl;
}
*/
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
    Eigen::Vector3d w1 = Eigen::Vector3d(0.0, 0.0, 1.0);
    Eigen::Vector3d v1 = Eigen::Vector3d(0.0, 0.0, 0.0);
    Eigen::Matrix4d e1 = math::matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Vector3d w2 = Eigen::Vector3d(0.0, 0.0, 1.0);
    Eigen::Vector3d v2 = Eigen::Vector3d(0.0, -L1, 0.0);
    Eigen::Matrix4d e2 = math::matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Vector3d w3 = Eigen::Vector3d(0.0, 0.0, 1.0);
    Eigen::Vector3d v3 = Eigen::Vector3d(0.0, -(L1+L2), 0.0);
    Eigen::Matrix4d e3 = math::matrix_exponential(w3, v3, joint_positions[2]);
    Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
    m.block<3, 1>(0, 3) = Eigen::Vector3d(L1+L2+L3, 0.0, 0.0);
    return e1 * e2 * e3 * m;
}

//-----------------------------------------------------------------------------------------------
// Task 5
// a)
/*Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions) // UR5
{
    //double L1 = 425.0, L2 = 392.0, H1 = 89.0, H2 = 95.0, W1 = 109.0, W2 = 82.0; // UR5
    Eigen::Vector3d w1 = Eigen::Vector3d(0.0, 0.0, 1.0);
    Eigen::Vector3d v1 = Eigen::Vector3d(0.0, 0.0, 0.0);
    Eigen::Matrix4d e1 = matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Vector3d w2 = Eigen::Vector3d(0.0, 1.0, 0.0);
    Eigen::Vector3d v2 = Eigen::Vector3d(-H1, 0.0, 0.0);
    Eigen::Matrix4d e2 = matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Vector3d w3 = Eigen::Vector3d(0.0, 1.0, 0.0);
    Eigen::Vector3d v3 = Eigen::Vector3d(-H1, 0.0, L1);
    Eigen::Matrix4d e3 = matrix_exponential(w3, v3, joint_positions[2]);
    Eigen::Vector3d w4 = Eigen::Vector3d(0.0, 1.0, 0.0);
    Eigen::Vector3d v4 = Eigen::Vector3d(-H1, 0.0, L1+L2);
    Eigen::Matrix4d e4 = matrix_exponential(w4, v4, joint_positions[3]);
    Eigen::Vector3d w5 = Eigen::Vector3d(0.0, 0.0, -1.0);
    Eigen::Vector3d v5 = Eigen::Vector3d(-W1, L1+L2, 0.0);
    Eigen::Matrix4d e5 = matrix_exponential(w5, v5, joint_positions[4]);
    Eigen::Vector3d w6 = Eigen::Vector3d(0.0, 1.0, 0.0);
    Eigen::Vector3d v6 = Eigen::Vector3d(H2-H1, 0.0, L1+L2);
    Eigen::Matrix4d e6 = matrix_exponential(w6, v6, joint_positions[5]);
    Eigen::Matrix4d m = Eigen::Matrix4d::Zero();
    m.block<3, 1>(0, 3) = Eigen::Vector3d(L1+L2, W1+W2, H1-H2);
    m.block<3, 3>(0, 0) << -1.0, 0.0, 0.0,
            0.0, 0.0, 1.0,
            0.0, 1.0, 0.0;
    m(3, 3) = 1.0;
    //print_pose("m", m);
    return e1 * e2 * e3 * e4 * e5 * e6 * m;
}*/
Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions)  //UR3e
{
    double L1 = 243.6, L2 = 213.2, H1 = 151.9, H2 = 85.4, W1 = 131.1, W2 = 92.1;
    Eigen::Vector3d w1 = Eigen::Vector3d(0.0, 0.0, 1.0);
    Eigen::Vector3d v1 = Eigen::Vector3d(0.0, 0.0, 0.0);
    Eigen::Matrix4d e1 = math::matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Vector3d w2 = Eigen::Vector3d(0.0, -1.0, 0.0);
    Eigen::Vector3d v2 = Eigen::Vector3d(H1, 0.0, 0.0);
    Eigen::Matrix4d e2 = math::matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Vector3d w3 = Eigen::Vector3d(0.0, -1.0, 0.0);
    Eigen::Vector3d v3 = Eigen::Vector3d(H1, 0.0, L1);
    Eigen::Matrix4d e3 = math::matrix_exponential(w3, v3, joint_positions[2]);
    Eigen::Vector3d w4 = Eigen::Vector3d(0.0, -1.0, 0.0);
    Eigen::Vector3d v4 = Eigen::Vector3d(H1, 0.0, L1+L2);
    Eigen::Matrix4d e4 = math::matrix_exponential(w4, v4, joint_positions[3]);
    Eigen::Vector3d w5 = Eigen::Vector3d(0.0, 0.0, -1.0);
    Eigen::Vector3d v5 = Eigen::Vector3d(W1, -L1-L2, 0.0);
    Eigen::Matrix4d e5 = math::matrix_exponential(w5, v5, joint_positions[4]);
    Eigen::Vector3d w6 = Eigen::Vector3d(0.0, -1.0, 0.0);
    Eigen::Vector3d v6 = Eigen::Vector3d(H1-H2, 0.0, L1+L2);
    Eigen::Matrix4d e6 = math::matrix_exponential(w6, v6, joint_positions[5]);
    Eigen::Matrix4d m = Eigen::Matrix4d::Zero();
    m.block<3, 1>(0, 3) = Eigen::Vector3d(-L1-L2, -W1-W2, H1-H2);
    m.block<3, 3>(0, 0) <<  1.0, 0.0, 0.0,
                                            0.0, 0.0, -1.0,
                                            0.0, 1.0, 0.0;
    m(3, 3) = 1.0;
    //print_pose("m", m);
    return e1 * e2 * e3 * e4 * e5 * e6 * m;
}

// b)
Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions) {
    double L1 = 243.6, L2 = 213.2, H1 = 151.9, H2 = 85.4, W1 = 131.1, W2 = 92.1;  // UR3e
    Eigen::Matrix4d t01;
    t01.block<3, 3>(0, 0) = math::rotate_z(joint_positions[0] * math::deg_to_rad);
    t01.block<3, 1>(0, 3) = Eigen::Vector3d(0.0, 0.0, H1);
    t01.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t01(3, 3) = 1.0;
    Eigen::Matrix4d t12;
    t12.block<3, 3>(0, 0) = math::rotate_y(-joint_positions[1] * math::deg_to_rad);
    t12.block<3, 1>(0, 3) = Eigen::Vector3d(0.0, -W1, 0.0);
    t12.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t12(3, 3) = 1.0;
    Eigen::Matrix4d t23;
    t23.block<3, 3>(0, 0) = math::rotate_y(-joint_positions[2] * math::deg_to_rad);
    t23.block<3, 1>(0, 3) = Eigen::Vector3d(-L1, 0.0, 0.0);
    t23.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t23(3, 3) = 1.0;
    Eigen::Matrix4d t34;
    t34.block<3, 3>(0, 0) = math::rotate_y(-joint_positions[3] * math::deg_to_rad);
    t34.block<3, 1>(0, 3) = Eigen::Vector3d(-L2, 0.0, 0.0);
    t34.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t34(3, 3) = 1.0;
    Eigen::Matrix4d t45;
    t45.block<3, 3>(0, 0) = math::rotate_z(-joint_positions[4] * math::deg_to_rad);
    t45.block<3, 1>(0, 3) = Eigen::Vector3d(0.0, 0.0, 0.0);
    t45.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t45(3, 3) = 1.0;
    Eigen::Matrix4d t56;
    t56.block<3, 3>(0, 0) = math::rotate_y(-joint_positions[5] * math::deg_to_rad);
    t56.block<3, 1>(0, 3) = Eigen::Vector3d(0.0, 0.0, -H2);
    t56.block<1, 3>(3, 0) = Eigen::Vector3d::Zero().transpose();
    t56(3, 3) = 1.0;
    Eigen::Matrix4d t67 = Eigen::Matrix4d::Identity();
    t67.block<3, 3>(0, 0) = math::rotate_x(90.0 * math::deg_to_rad);
    t67.block<3, 1>(0, 3) = Eigen::Vector3d(0.0, -W2, 0.0);
    Eigen::Matrix4d m = Eigen::Matrix4d::Zero();
    m.block<3, 1>(0, 3) = Eigen::Vector3d(-L1 - L2, -W1 - W2, H1 - H2);
    m.block<3, 3>(0, 0) <<  1.0, 0.0, 0.0,
                                            0.0, 0.0, -1.0,
                                            0.0, 1.0, 0.0;
    m(3, 3) = 1.0;
    //print_pose("m", m);
    return t01 * t12 * t23 * t34 * t45 * t56 * t67;
}





int main()
{
    /*Eigen::Vector3d e = Eigen::Vector3d(60.0, 45.0, 30.0) * math::deg_to_rad;
    Eigen::Matrix3d rotation_matrix= math::rotation_matrix_from_euler("zyx", e);
    Eigen::Vector3d ea = euler_zyx_from_rotation(rotation_matrix);
    std::cout << " E:" << e.transpose() * math::rad_to_deg << std::endl;
    std::cout << "Ea:" << ea.transpose() * math::rad_to_deg << std::endl;
    task2_a();
    task2_b();
    print_pose("shit", Eigen::Matrix4d::Identity());
     */
    //task4b-c
    std::vector<double> joint_positions = {10.0, -15.0, 2.75};
    Eigen::Matrix4d transform = planar_3r_fk_transform(joint_positions);
    Eigen::Matrix4d screw = planar_3r_fk_screw(joint_positions);
    math::print_pose("end position transform:", transform);
    math::print_pose("end position screw", screw);
    //task5a
    //std::vector<double> joint_positions = {30.0, -60.0, 30.0, -50.0, 90.0, 0.0};  //-286, -317, 545, & 9.9, 0, -59,8
    std::vector<double> joint_positions1 = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  //-457,  -223,  67 & 90, 0, 0
    //std::vector<double> joint_positions = {0.0, 0.0, 0.0, -90.0, 0.0, 0.0};  //-542, -223, 152 & -1.3, 89.8, -91
    //std::vector<double> joint_positions = {0.0, 180.0, 0.0, 0.0, 0.0, 0.0};  //456, -223, 236 & -90, 0, 180
    //std::vector<double> joint_positions = {0.0, 90, 0.0, 0.0, 0.0, 0.0};  //-85, -223, 608 & -178.2, 90, 92
    //std::vector<double> joint_positions = {0.0, 0.0, 0.0, 0.0, 180.0, 0.0};  //-457, -40, 67 & 90, 0, -180
    //std::vector<double> joint_positions1 = {0.0, 0.0, 0.0, 90.0, 0.0, 0.0};  //-371, -223, 152 & 2, -90, 88
    //std::vector<double> joint_positions = {0.0, 0.0, 0.0, 30.0, 0.0, 0.0};  //-414, -223, 79 & 90, -30, 0
    //std::vector<double> joint_positions = {0.0, 0.0, 0.0, -30.0, 0.0, 0.0};  //-499, -223, 79 & 90, 30, 0


    Eigen::Matrix4d screw_ur3e = ur3e_fk_screw(joint_positions1);
    Eigen::Matrix4d transform_ur3e = ur3e_fk_transform(joint_positions1);
    math::print_pose("end position screw ur3e:", screw_ur3e);
    math::print_pose("end position transform ur3e:", transform_ur3e);


    /*for (int i = 0; i < 180; ++i) {
     joint_positions[4] = -90 + i;
     Eigen::Matrix4d screw_ur3e = ur3e_fk_screw(joint_positions);

     // Extract the 3x3 rotation matrix
     Eigen::Matrix3d R = screw_ur3e.block<3, 3>(0, 0);
     Eigen::Vector3d e = euler_zyx_from_rotation(R)*rad_to_deg;

     std::cout << "i = " << i << "      " << e.transpose() << std::endl;
     //std::cout << e.transpose() << std::endl;

 }*/
    return 0;
}