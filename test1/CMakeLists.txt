find_package(Eigen3 CONFIG REQUIRED)

add_executable(test1 main.cpp)

target_link_libraries(test1 PRIVATE Eigen3::Eigen)
