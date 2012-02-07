cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)

project(ParameterizedCFM) 

# ITK 4.0 required
FIND_PACKAGE(ITK 4.0 REQUIRED) 
INCLUDE(${ITK_USE_FILE})

SET (EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin CACHE INTERNAL 
"Single output directory for building all executables.")

INCLUDE_DIRECTORIES(
  ${CMAKE_SOURCE_DIR}/Utilities/vul
  ${CMAKE_SOURCE_DIR}/Utilities/json
  )

SUBDIRS(Utilities Code)
