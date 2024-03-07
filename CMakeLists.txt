cmake_minimum_required(VERSION 3.16)
project(example)
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Libigl
include(libigl)
find_package(CGAL 5.5 EXACT REQUIRED PATHS /home/azhar/Downloads)
find_package(Boost)
find_package(GMP REQUIRED)

# Enable the target igl::glfw
igl_include(glfw)

# Add your project files
file(GLOB SRC_FILES *.cpp)
add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC igl::glfw ${CGAL_LIBRARIES} ${GMP_LIBRARIES}) 

