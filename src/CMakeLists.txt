cmake_minimum_required(VERSION 3.0.0)
project(intel-spmxv VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})
link_directories(${MKL_LIBRARIES})

add_executable(intel-spmxv intel-spmxv.c)

add_subdirectory(deps)

target_sources(intel-spmxv PRIVATE
    mm-utils.c
)

target_include_directories(intel-spmxv PRIVATE
    include
)

find_package(OpenMP REQUIRED)

target_link_libraries(intel-spmxv PRIVATE
    mkl_intel_lp64
    mkl_sequential
    mkl_core
    OpenMP::OpenMP_C
    mmio
    m
)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
