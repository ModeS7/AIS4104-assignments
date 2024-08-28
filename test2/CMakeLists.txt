find_package(Eigen3 CONFIG REQUIRED)

add_executable(test2 main.cpp)

target_link_libraries(test2 PRIVATE Eigen3::Eigen)
