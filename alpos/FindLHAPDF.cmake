# - Try to find LHAPDF
# Once done, this will define
#
# LHAPDF_FOUND
# LHAPDF_INCLUDE_DIRS
# LHAPDF_LIBRARIES


MESSAGE(STATUS "Looking for LHAPDF...")

# try to use lhapdf-config tool
FIND_PROGRAM(LHAPDF_CONFIG_EXECUTABLE NAMES lhapdf-config)
EXEC_PROGRAM(${LHAPDF_CONFIG_EXECUTABLE} ARGS "--prefix" OUTPUT_VARIABLE LHAPDF_CONFIG__PREFIX)

FIND_LIBRARY(LHAPDF_LIBRARIES
  NAMES LHAPDF
  PATHS
  ${LHAPDF_PREFIX}/lib
  ${LHAPDF_CONFIG_PREFIX}/lib
  ${PREFIX}/lib
  $ENV{HOME}/.local/lib
  /usr/local/lib
  /usr/lib
  )

FIND_PATH(LHAPDF_INCLUDE_DIRS
  LHAPDF/LHAPDF.h
  $ENV{HOME}/.local/include
  /usr/include
  /usr/local/include
  ${LHAPDF_PREFIX}/include
  ${LHAPDF_CONFIG_PREFIX}/include
  )

if (NOT LHAPDF_LIBRARIES)
  MESSAGE(STATUS "LHAPDF not found.")
  MESSAGE(STATUS "You can pass -DLHAPDF_PREFIX=/path/to/lhapdf to cmake.")
else()
  MESSAGE(STATUS "LHAPDF found. Library path: " ${LHAPDF_LIBRARIES}  ", Include path:" ${LHAPDF_INCLUDE_DIRS})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LHAPDF DEFAULT_MSG LHAPDF_LIBRARIES LHAPDF_INCLUDE_DIRS)
