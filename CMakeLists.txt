cmake_minimum_required(VERSION 3.16)
project(gpod)

math(EXPR stack_size "100*1024*1024")
set(CMAKE_EXE_LINKER_FLAGS "-Wl,--stack,${stack_size}")
set(CMAKE_CXX_FLAGS_RELEASE  "-O2")

add_subdirectory(lib/sofa)
add_subdirectory(lib/nrlmsise)
add_subdirectory(lib/openblas)
add_subdirectory(lib/gsl)
add_subdirectory(lib/boost)

include_directories(${PROJECT_BINARY_DIR})
include_directories(src)
include_directories(lib)
include_directories(lib/sofa/src)
include_directories(lib/nrlmsise)

link_directories(lib)

aux_source_directory(src SRC_GPOD)

if(WIN32)
    set(LIBS Winmm)
else()
    set(LIBS m pthread)
endif()

add_library(gpod ${SRC_GPOD})
target_link_libraries(gpod ${LIBS} m openblas gsl sofa nrlmsise pthread odeint)

add_executable(test_podflt test_podflt.c)
target_link_libraries(test_podflt gpod)








