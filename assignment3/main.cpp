#include <iostream>

#include <Eigen/Dense>
#include "math/math.h"
#include <fstream>


using namespace math;


// Task 1
// a)
Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v)
{
    Eigen::VectorXd eigen_vector(v.size());
    for (int i = 0; i < v.size(); i++)
    {
        eigen_vector(i) = v[i];
    }
    return eigen_vector;
}

// b)
bool is_average_below_eps(const std::vector<double> &values, double eps = 10e-7, uint8_t n_values = 5u)
{
    double sum = 0.0;
    for (int i = 0; i < values.size(); i++)
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
{
    double L1 = 243.6, L2 = 213.2, H1 = 151.9, H2 = 85.4, W1 = 131.1, W2 = 92.1;
    std::vector<Eigen::Vector3d> s{
        Eigen::Vector3d(0.0, 0.0, 1.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, -1.0),
        Eigen::Vector3d(0.0, 1.0, 0.0)
    };
    std::vector<Eigen::Vector3d> q{
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, H1),
        Eigen::Vector3d(L1, 0.0, H1),
        Eigen::Vector3d(L1+L2, 0.0, H1),
        Eigen::Vector3d(L1+L2, W1, 0.0),
        Eigen::Vector3d(L1+L2, 0.0, H1-H2)
    };

    std::vector<Eigen::VectorXd> screws(s.size());
    for (int i = 0; i < s.size(); i++)
    {
        screws[i] = math::screw_axis(q[i], s[i], 0.0);
    }

    Eigen::Matrix4d m = Eigen::Matrix4d::Identity();
    m.block<3, 1>(0, 3) = Eigen::Vector3d(L1+L2, W1+W2, H1-H2);
    m.block<3, 3>(0, 0) << -1.0, 0.0, 0.0,
                                            0.0, 0.0, 1.0,
                                            0.0, 1.0, 0.0;
    return std::make_pair(m, screws);
}

// d)
Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd &screw, double theta)
{
    return math::matrix_exponential(screw.head(3), screw.tail(3), theta);
}

Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions)
{
    /*
     * Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
     * sdt::vector<Eigen::VectorXd> screws;
     * std::tie(M, screws) = ur3e_space_chain();
     */
    auto [m, space_screws] = ur3e_space_chain();
    Eigen::Matrix4d t06 = Eigen::Matrix4d::Identity();
    for (int i = 0; i < joint_positions.size(); i++){
        t06 *= matrix_exponential(space_screws[i], joint_positions[i]);
    }
    return t06 * m;
}

// e)
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain()
{
    auto [m_sb, space_screws] = ur3e_space_chain();
    Eigen::Matrix4d m_bs = m_sb.inverse();
    std::vector<Eigen::VectorXd> body_twists(space_screws.size());
    Eigen::MatrixXd adj = math::adjoint_matrix(m_bs);
    for (int i = 0; i < space_screws.size(); i++)
    {
        body_twists[i] = adj * space_screws[i];
        //std::cout << "space_screws: " << space_screws[i] << std::endl;
        //std::cout << "screw: " << screws[i] << std::endl;
    }
    return std::make_pair(m_sb, body_twists);
}

// f)
Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions)
{
    auto [m_sb, body_twist] = ur3e_body_chain();
    Eigen::Matrix4d t06 = m_sb;
    for (int i = 0; i < joint_positions.size(); i++){
        t06 *= matrix_exponential(body_twist[i], joint_positions[i]);
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
                                                     double x_0, double dx_0 = 0.5, double eps = 10e-7)
{
    double x = x_0;
    double dx = dx_0;
    int iterations = 0;
    while (std::abs(f(x)) > eps)
    {
        x = x - f(x)/((f(x+dx)-f(x))/dx);
        iterations++;
        if (iterations > 10000){
            std::cout << "reatched max iterations of:" << iterations++ << std::endl;
            break;
        }
    }
    return std::make_pair(iterations, x);
}
/*void newton_r_r_f_test() {
    std::cout << "newton_raphson_root_find tests" << std::endl;
    auto [iter, x] = newton_raphson_root_find([](double x) { return x * x - 4.0; }, 1.0);
    std::cout << "Iterations: " << iter << " Root: " << x << std::endl;
}*/

// b)
std::pair<uint32_t, double> gradient_descent_root_find(const std::function<double(double)> &f,
                                                       double x_0, double gamma = 0.1, double dx_0 = 0.5, double eps = 1e-7)
{
    double x_old = x_0;
    double x = (x_0 - gamma * ((f(x_0 + dx_0) - f(x_0)) / dx_0));
    double error_new = {std::abs(f(x))};
    double error = {std::abs(f(x_old))};
    double x_new;
    int iterations = 1;

    while (error_new > eps){
        x_new = (x - gamma * (error_new - error) / (x - x_old));
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
{
    auto [m, space_screws] = ur3e_space_chain();
    Eigen::MatrixXd jacobian(6, space_screws.size());
    int n_c = current_joint_positions.size();
    int n_s = space_screws.size();
    if (n_c != n_s){
        std::cerr << "Invalid number of joint positions" << std::endl;
        jacobian = Eigen::MatrixXd::Zero(6, n_s);
    }
    else{
        Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
        for (int i = 0; i < n_s; i++)
        {
            if (i == 0){
                jacobian.block<6, 1>(0, i) = space_screws[i];
                continue;
            }
            else{
                Eigen::MatrixXd adj;
                t *= matrix_exponential(space_screws[i-1], current_joint_positions[i-1]);
                adj = math::adjoint_matrix(t);
                jacobian.block<6, 1>(0, i) = adj * space_screws[i];
            }
        }
        return jacobian;
    }
}

// b)
Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    auto [m, body_twists] = ur3e_body_chain();
    Eigen::MatrixXd jacobian(6, body_twists.size());
    int n_c = current_joint_positions.size();
    int n_s = body_twists.size();
    if (n_c != n_s){
        std::cerr << "Invalid number of joint positions" << std::endl;
        jacobian = Eigen::MatrixXd::Zero(6, n_s);
    }
    else{
        Eigen::Matrix4d t = Eigen::Matrix4d::Identity();
        for (int i = 0; i < n_s; i++)
        {
            Eigen::MatrixXd adj;
            t *= matrix_exponential(body_twists[i], current_joint_positions[i]);
            adj = math::adjoint_matrix(t);
            jacobian.block<6, 1>(0, i) = adj * body_twists[i];
        }
        return jacobian;
    }
}



int main()
{
    // Task 1
    //ur3e_test_fk();
    // Task 2
    test_root_find();
    return 0;
}
























