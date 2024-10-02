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
    Eigen::Matrix3d mr = math::rotate_y(-90.0 * math::deg_to_rad)
            * math::rotate_x(-90.0 * math::deg_to_rad)
            * math::rotate_z(-90.0 * math::deg_to_rad);
    Eigen::Matrix4d m = math::transformation_matrix(mr, Eigen::Vector3d{L1 + L2, W1 + W2, H1 - H2});

    std::vector<Eigen::VectorXd> screws{
        math::screw_axis({0.0, 0.0, 0.0},{0.0, 0.0, 1.0}, 0.0),
        math::screw_axis({0.0, 0.0, H1},{0.0, 1.0, 0.0}, 0.0),
        math::screw_axis({L1, 0.0, H1},{0.0, 1.0, 0.0}, 0.0),
        math::screw_axis({L1+L2, 0.0, H1},{0.0, 1.0, 0.0}, 0.0),
        math::screw_axis({L1+L2, W1, 0.0},{0.0, 0.0, -1.0}, 0.0),
        math::screw_axis({L1+L2, 0.0, H1-H2},{0.0, 1.0, 0.0}, 0.0)
    };
    return std::make_pair(m, screws);
}



int main()
{

    return 0;
}