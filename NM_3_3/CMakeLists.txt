cmake_minimum_required(VERSION 3.24)
project(NM_3_3)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES main.cpp)
set(PYTHON_FILES plot.py)

add_executable(NM_3_3 ${SOURCE_FILES} ${PYTHON_FILES})
target_link_libraries(NM_3_3 ${PYTHON_LIBRARIES})
