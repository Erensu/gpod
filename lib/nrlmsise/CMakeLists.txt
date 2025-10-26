
set(CMAKE_CXX_FLAGS_RELEASE  "-O2")

include_directories(.)

add_library(nrlmsise eopspw.c nrlmsise-00.c nrlmsise-00_data.c)
target_link_libraries(nrlmsise m)

