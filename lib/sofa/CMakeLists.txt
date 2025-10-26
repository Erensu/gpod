
set(CMAKE_CXX_FLAGS_RELEASE  "-O2")

aux_source_directory(src SRC_SOFA)
include_directories(src)

add_library(sofa ${SRC_SOFA})
target_link_libraries(sofa m)

