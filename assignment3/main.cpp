#include <iostream>

#include <Eigen/Dense>
#include "math/math.h"
#include <fstream>


using namespace math;


// Task 1
// a)
Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v)
// Convert a std::vector to an Eigen::VectorXd
{
    Eigen::VectorXd eigen_vector(v.size());
    for (size_t i = 0; i < v.size(); i++)
        eigen_vector(i) = v[i];
    return eigen_vector;
}

// b)
bool is_average_below_eps(const std::vector<double> &values, double eps = 10e-7, uint8_t n_values = 5u)
// Check if the average of the values is below eps
{
    double sum = 0.0;
    for (size_t i = 0; i < values.size(); i++)
    {
        sum += values[i];
    }
    double average = sum / values.size();
    if (values.size() < n_values)
    {
        return false;
    }
    else if (average < eps)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// c)
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_space_chain()
// Return the space chain of the UR3e robot
// Modern Robotics page 146.(UR5)
{
    const double L1 = 0.2436, L2 = 0.2132, H1 = 0.1519, H2 = 0.0854, W1 = 0.1311, W2 = 0.0921;
    const std::vector<Eigen::Vector3d> s{
        {0.0, 0.0, 1.0},
        {0.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, -1.0},
        {0.0, 1.0, 0.0}
    };
    const std::vector<Eigen::Vector3d> q{
        {0.0, 0.0, 0.0},
        {0.0, 0.0, H1},
        {L1, 0.0, H1},
        {L1+L2, 0.0, H1},
        {L1+L2, W1, 0.0},
        {L1+L2, 0.0, H1-H2}
    };

    std::vector<Eigen::VectorXd> screws(s.size());
    for (size_t i = 0; i < s.size(); i++)
    {
        screws[i] = math::screw_axis(q[i], s[i], 0.0);
    }

    Eigen::Matrix4d m;
    m <<
    1.0, 0.0, 0.0, L1+L2,
    0.0, 0.0, 1.0, W1+W2,
    0.0, -1.0, 0.0, H1-H2,
    0.0, 0.0, 0.0, 1.0;
    return std::make_pair(m, screws);
}

// d)
Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions)
// Compute the forward kinematics of the UR3e robot in space frame.
// Modern Robotics page 146.
{
    Eigen::Matrix4d m;
    Eigen::Matrix4d t06 = Eigen::Matrix4d::Identity();
    std::vector<Eigen::VectorXd> space_screws;
    std::tie(m, space_screws) = ur3e_space_chain();
    for (size_t i = 0; i < joint_positions.size(); i++){
        t06 *= math::matrix_exponential(space_screws[i].head(3), space_screws[i].tail(3), joint_positions[i]);
    }
    return t06 * m;
}

// e)
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain()
// Return the body chain of the UR3e robot converting
// the space chain to the body chain. Modern Robotics page 99
{
    Eigen::Matrix4d m_sb;
    std::vector<Eigen::VectorXd> space_screws;
    std::tie(m_sb, space_screws) = ur3e_space_chain();
    const Eigen::Matrix4d m_bs = m_sb.inverse();
    std::vector<Eigen::VectorXd> body_twists(space_screws.size());
    Eigen::MatrixXd adj = math::adjoint_matrix(m_bs);
    for (size_t i = 0; i < space_screws.size(); i++)
    {
        body_twists[i] = adj * space_screws[i];
    }
    return std::make_pair(m_sb, body_twists);
}

// f)
Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions)
// Compute the forward kinematics of the UR3e robot in body frame.
// Modern Robotics page 148.
{
    Eigen::Matrix4d m_sb;
    std::vector<Eigen::VectorXd> body_twist;
    std::tie(m_sb, body_twist) = ur3e_body_chain();
    Eigen::Matrix4d t06 = m_sb;
    for (size_t i = 0; i < joint_positions.size(); i++){
        t06 *= math::matrix_exponential(body_twist[i].head(3), body_twist[i].tail(3), joint_positions[i]);
    }
    return t06;
}

// g)
void ur3e_test_fk()
{
    std::cout << "Forward kinematics tests" << std::endl;
    math::print_pose("pos1_s", ur3e_space_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0})));
    math::print_pose("pos1_b", ur3e_body_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0})));
    std::cout << std::endl;
    math::print_pose("pos2_s", ur3e_space_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0})));
    math::print_pose("pos2_b", ur3e_body_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0})));
    std::cout << std::endl;
    math::print_pose("pos3_s", ur3e_space_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0})));
    math::print_pose("pos3_b", ur3e_body_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0})));
    std::cout << std::endl;
    math::print_pose("pos4_s", ur3e_space_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0})));
    math::print_pose("pos4_b", ur3e_body_fk(std_vector_to_eigen(
            std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0})));
}

// Task 2
// a)
std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f,
                                                     double x_0,
                                                     double dx_0 = 0.5,
                                                     double eps = 10e-7)
// Find the root of f using the Newton-Raphson method. Modern Robotics page 225.
{
    uint32_t max_iterations = 10000;
    double x = x_0;
    double dx = dx_0;
    uint32_t iterations = 0;
    double fx = f(x);
    while (std::abs(f(x)) > eps)
    {
        double dfdx = (f(x+dx)-fx)/dx;
        if (std::abs(dfdx) < eps){
            std::cout << "Derivative is zero, reatched iteration: " << iterations++ << std::endl;
            break;
        }
        x = x - fx/dfdx;
        fx = f(x);
        iterations++;
        if (iterations >= max_iterations){
            std::cout << "reatched max iterations of:" << iterations++ << std::endl;
            break;
        }
    }
    return std::make_pair(iterations, x);
}

// b)
std::pair<uint32_t, double> gradient_descent_root_find(const std::function<double(double)> &f,
                                                       double x_0,
                                                       double gamma = 0.01,
                                                       double dx_0 = 0.5,
                                                       double eps = 1e-7)
// Find the root of f using the gradient descent method minimizing error.
{
    double x_old = x_0;
    double x = (x_0 - gamma * ((f(x_0 + dx_0) - f(x_0)) / dx_0));
    double error_new = {std::abs(f(x))};
    double error = {std::abs(f(x_old))};
    double x_new;
    uint32_t iterations = 1;

    while (error_new > eps){
        double dedx = (error_new - error) / (x - x_old);
        if (std::abs(dedx) < eps){
            std::cout << "Derivative is zero, reatched iteration: " << iterations++ << std::endl;
            break;
        }
        x_new = (x - gamma * dedx);
        error = error_new;
        error_new = (std::abs(f(x_new)));
        x_old = x;
        x = x_new;
        iterations++;

        if (iterations >= 100000) {
            std::cout << "Reached maximum iterations: " << iterations << std::endl;
            break;
        }
    }
    return std::make_pair(iterations, x);
}

// c)
void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0)
{
    auto [iterations, x_hat] = newton_raphson_root_find(f, x0);
    std::cout << "NR root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) << std::endl;
}
void test_gradient_descent_root_find(const std::function<double(double)> &f, double x0)
{
    auto [iterations, x_hat] = gradient_descent_root_find(f, x0);
    std::cout << "GD root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" << f(x_hat) << std::endl;
}
void test_root_find()
{
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x)
    {
        return (x - 3.0) * (x - 3.0) - 1.0;
    };
    test_newton_raphson_root_find(f1, -20);
    test_gradient_descent_root_find(f1, -20);
}

// Task 3
// a)
Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions)
// Compute the space jacobian of the UR3e robot. Modern Robotics page 178.
{
    Eigen::Matrix4d m;
    std::vector<Eigen::VectorXd> space_screws;
    std::tie(m, space_screws) = ur3e_space_chain();
    const int n_c = current_joint_positions.size();
    const int n_s = space_screws.size();
    if (n_c != n_s){
        std::cerr << "Invalid number of joint positions" << std::endl;
        return Eigen::MatrixXd::Zero(6, n_s);;
    }
    Eigen::MatrixXd jacobian(6, space_screws.size());
    Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
    size_t i = 0;
    jacobian.block<6, 1>(0, i) = space_screws[i];
    i++;
    while (i < n_s){
        Eigen::MatrixXd adj;
        t *= math::matrix_exponential(space_screws[i-1].head(3), space_screws[i-1].tail(3), current_joint_positions[i-1]);
        adj = math::adjoint_matrix(t);
        jacobian.block<6, 1>(0, i) = adj * space_screws[i];
        i++;
    }
    return jacobian;
}

// b)
Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions)
// Compute the body jacobian of the UR3e robot. Modern Robotics page 181.
{
    Eigen::Matrix4d m;
    std::vector<Eigen::VectorXd> body_twists;
    std::tie(m, body_twists) = ur3e_body_chain();
    Eigen::MatrixXd jacobian(6, body_twists.size());
    const int n_c = current_joint_positions.size();
    const int n_s = body_twists.size();
    Eigen::MatrixXd adj;
    Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
    if (n_c != n_s){
        std::cerr << "Invalid number of joint positions" << std::endl;
        return Eigen::MatrixXd::Zero(6, n_s);
    }
    size_t i = n_s-1;
    jacobian.block<6, 1>(0, i) = body_twists[i];
    i--;
    while (i+1 > 0)
    {
        Eigen::VectorXd b_t = -body_twists[i+1];
        t *= matrix_exponential(b_t.head(3), b_t.tail(3), current_joint_positions[i+1]);
        adj = math::adjoint_matrix(t);
        jacobian.block<6, 1>(0, i) = adj * body_twists[i];
        i--;
    }
    return jacobian;

}

// c)
void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions)
{
    Eigen::Matrix4d tsb = ur3e_body_fk(joint_positions);
    auto [m, space_screws] = ur3e_space_chain();
    Eigen::MatrixXd jb = ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = adjoint_matrix(tsb.inverse());
    std::cout << "Jb: " << std::endl << jb << std::endl << "Ad_tbs*Js:" << std::endl << ad_tbs * js << std::endl << std::endl;
    std::cout << "Js: " << std::endl << js << std::endl << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb << std::endl << std::endl;
    std::cout << "d Jb: " << std::endl << jb - ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js - ad_tsb * jb << std::endl << std::endl;
}
void ur3e_test_jacobian()
{
    std::cout << "Jacobian matrix tests" << std::endl;
    ur3e_test_jacobian(std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}));
    ur3e_test_jacobian(std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}));
}

// Task 4
// a)
std::pair<size_t, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd,
                                                const Eigen::VectorXd &current_joint_positions,
                                                double gamma = 1e-2,
                                                double v_e = 4e-3,
                                                double w_e = 4e-3)
// Compute the joint positions of the UR3e robot
// using the body jacobian. Modern Robotics page 228-229.
{
    Eigen::Matrix4d t_sb = ur3e_space_fk(current_joint_positions);
    Eigen::Matrix4d t_err = (t_sb).inverse() * t_sd;
    Eigen::VectorXd v_b;
    double theta;
    std::tie(v_b, theta) = math::matrix_logarithm(t_err);
    v_b *= theta;
    size_t iterations = 0;
    Eigen::VectorXd joint_p = current_joint_positions;
    Eigen::MatrixXd jacobian_b;
    while ((v_b.head(3).norm() > v_e) || (v_b.tail(3).norm() > w_e))
    {
        jacobian_b = ur3e_body_jacobian(joint_p);
        joint_p += gamma * jacobian_b.completeOrthogonalDecomposition().pseudoInverse() * v_b;
        t_sb = ur3e_space_fk(joint_p);
        t_err = t_sb.inverse() * t_sd;
        std::tie(v_b, theta) = math::matrix_logarithm(t_err);
        v_b *= theta;
        iterations++;
        if (iterations >= 10000){
            std::cout << "Reached max iterations: " << iterations << std::endl;
            break;
        }
    }
    return std::make_pair(iterations, joint_p);
}

// b)
void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd &j0)
{
std::cout << "Test from pose" << std::endl;
Eigen::Matrix4d t_sd = transformation_matrix(rotation_matrix_from_euler("zyx", zyx), pos);
auto [iterations, j_ik] = ur3e_ik_body(t_sd, j0);
Eigen::Matrix4d t_ik = ur3e_body_fk(j_ik);
print_pose(" IK pose",t_ik);
print_pose("Desired pose", t_sd);
std::cout << "Converged after " << iterations << " iterations" << std::endl;
std::cout << "J_0: " << j0.transpose() * math::rad_to_deg << std::endl;
std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg << std::endl << std::endl;
}
void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0)
{
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = ur3e_space_fk(joint_positions);
    auto [iterations, j_ik] = ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = ur3e_body_fk(j_ik);
    print_pose(" IK pose", t_ik);
    print_pose("Desired pose", t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() * math::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * math::rad_to_deg << std::endl << std::endl;
}
void ur3e_ik_test()
{
    Eigen::VectorXd j_t0 = std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) * math::deg_to_rad;
    Eigen::VectorXd j_t1 = std_vector_to_eigen(std::vector<double>{0.0, 0.0, -89.0, 0.0, 0.0, 0.0}) * math::deg_to_rad;
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} * math::deg_to_rad, j_t0);
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} * math::deg_to_rad, j_t1);
    Eigen::VectorXd j_t2 = std_vector_to_eigen(std::vector<double>{50.0, -30.0, 20, 0.0, -30.0, 50.0}) * math::deg_to_rad;
    Eigen::VectorXd j_d1 = std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) * math::deg_to_rad;
    ur3e_ik_test_configuration(j_d1, j_t0);
    ur3e_ik_test_configuration(j_d1, j_t2);
}


int main()
{
    // Task 1
    ur3e_test_fk();
    // Task 2
    test_root_find();
    // Task 3
    ur3e_test_jacobian();
    // Task 4
    ur3e_ik_test();
    return 0;
}

