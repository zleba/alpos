# - Try to find fastNLO
# Once done, this will define
#
# FNLO_FOUND
# FNLO_INCLUDE_DIRS
# FNLO_LIBRARIES

MESSAGE(STATUS "Find fastNLO...")

FIND_LIBRARY(FNLO_LIBRARIES
  NAMES fastnlotoolkit
  PATHS
  ${FNLO_PREFIX}/lib
  ${PREFIX}/lib
  $ENV{HOME}/.local/lib
  /usr/local/lib
  /usr/lib
  )

FIND_PATH(FNLO_INCLUDE_DIRS
  NAMES fastnlotk/fastNLOReader.h
  PATHS
  ${FNLO_PREFIX}/include
  $ENV{HOME}/.local/include
  /usr/local/include
  /usr/include
  )

if (NOT FNLO_LIBRARIES)
  MESSAGE(STATUS "FNLO not found.")
  MESSAGE(STATUS "You can pass -DFNLO_PREFIX=/path/to/lhapdf to cmake.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FastNLO DEFAULT_MSG FNLO_LIBRARIES FNLO_INCLUDE_DIRS)
