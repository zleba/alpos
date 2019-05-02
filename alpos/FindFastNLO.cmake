# - Try to find fastNLO
# Once done, this will define
#
# FASTNLO_FOUND
# FASTNLO_INCLUDE_DIRS
# FASTNLO_LIBRARIES

MESSAGE(STATUS "Find fastNLO...")

FIND_LIBRARY(FASTNLO_LIBRARIES
  NAMES fastnlotoolkit
  PATHS
  ${FASTNLO_PREFIX}/lib
  ${PREFIX}/lib
  $ENV{HOME}/.local/lib
  /usr/local/lib
  /usr/lib
  )

FIND_PATH(FASTNLO_INCLUDE_DIRS
  NAMES fastnlotk/fastNLOReader.h
  PATHS
  ${FASTNLO_PREFIX}/include
  $ENV{HOME}/.local/include
  /usr/local/include
  /usr/include
  )

if (NOT FASTNLO_LIBRARIES)
  MESSAGE(STATUS "FASTNLO not found.")
  MESSAGE(STATUS "You can pass -DFASTNLO_PREFIX=/path/to/lhapdf to cmake.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FastNLO DEFAULT_MSG FASTNLO_LIBRARIES FASTNLO_INCLUDE_DIRS)
