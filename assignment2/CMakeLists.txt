find_package(Eigen3 CONFIG REQUIRED)

add_executable(assignment2 main.cpp)

target_link_libraries(assignment2 PUBLIC math PRIVATE Eigen3::Eigen)
