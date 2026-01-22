cmake_minimum_required(VERSION 3.5)

include_directories(htslib/include)
link_directories(htslib/lib)

project(anomalous_junction_count6 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(anomalous_junction_count6 main.cpp anomalous_junction_count.cpp anomalous_junction_count.h gtf_file.cpp gtf_file.h bam_file.cpp bam_file.h interval_tree.h progress_bar.h)

find_package(Threads REQUIRED)
target_link_libraries (anomalous_junction_count6 PRIVATE libhts.a  z m bz2 lzma curl crypto deflate Threads::Threads)
