cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

project(PNP_FullResidual_MMSGF)

add_executable(pnp
  # header files (not needed, but added so they show up in IDEs)
        include/PNPGridField.h
        include/PNPInputData.h
        include/PNPNodeData.h
        include/PNPManufacturedSoln.h
        include/UtilFunctions.h
        include/PNPBaseEquation.h
  # source files to compile
  src/main.cpp
)

# add the ./include directory to the include search path so we can #include
# the header files we put in include/
target_include_directories(pnp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

if (NOT TARGET talyfem)
  find_package(talyfem REQUIRED)
endif()

# this will also copy all the compiler definitions, etc. that TalyFEM was
# compiled with (that were marked as PUBLIC) - e.g. this will automatically
# add "-DENABLE_HDF5" if TalyFEM was compiled with HDF5 support
target_link_libraries(pnp talyfem)
