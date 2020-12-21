cmake_minimum_required(VERSION 3.1)
project(fluid)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
include_directories("include/")

file(GLOB SRCFILES src/*.cpp)
add_library(core ${SRCFILES})
target_link_libraries(core igl::core igl::opengl igl::opengl_glfw )

add_executable(fluid "main.cpp")
target_link_libraries(${PROJECT_NAME} core igl::core igl::opengl_glfw)
#target_link_libraries(example core igl::core igl::opengl igl::opengl_glfw )


