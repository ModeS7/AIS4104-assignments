#include <iostream>

#include <Eigen/Dense>

double deg_to_rad(double degrees)
{
    return degrees * 0.0174533;
}

double rad_to_deg(double radians)
{
    return radians * 57.2957795;
}

Eigen::Matrix3d rotate_x(double angle)
{
    Eigen::Matrix3d matrix;{
             1.0, 0.0, 0.0,
            0.0, std::cos(angle), std::sin(angle),
            0.0, std::sin(angle), std::cos(angle);
    }
    return matrix;
}



void example(double constant)
{
    Eigen::Matrix3d identity;
    identity <<
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
    std::cout << "I: " << std::endl << identity << std::endl << std::endl;
    std::cout << constant <<"*I: " << std::endl << constant * identity << std::endl << std::endl;
}

int main()
{
    example(2.0);
    Eigen::Matrix3d x_0 = rotate_x(angle:0.0);
    std::cout <<x_0 << std::endl;

    return 0;
}
