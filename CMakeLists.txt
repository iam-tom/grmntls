cmake_minimum_required (VERSION 2.6)
project(Garmin)
include_directories(include /opt/local/include /home/tom/workspace_cpp/eigen/)
link_directories (${DIR}/src/utils /opt/local/lib ${DIR}/lib)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(EXECUTABLE_OUTPUT_PATH bin)
set(LIBRARY_OUTPUT_PATH lib)


FIND_PACKAGE(Boost COMPONENTS system filesystem REQUIRED)
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR}  )

add_library(usb_comm src/usb_comm.c)

add_executable(mytest  src/tst.cpp)
target_link_libraries(mytest usb_comm garmintools usb ${Boost_LIBRARIES})

set(CMAKE_BUILD_TYPE Debug)
