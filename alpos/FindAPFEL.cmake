# - Try to find Apfel
# Once done, this will define
#
# APFEL_FOUND
# APFEL_INCLUDE_DIRS
# APFEL_LIBRARIES

MESSAGE(STATUS "Find Apfel...")

FIND_LIBRARY(APFEL_LIBRARIES
  NAMES APFEL
  PATHS
  ${APFEL_PREFIX}/lib
  ${PREFIX}/lib
  /usr/local/lib
  /usr/lib
  )

FIND_PATH(APFEL_INCLUDE_DIRS
  APFEL/APFEL.h
  /usr/include
  /usr/local/include
  ${APFEL_PREFIX}/include
  )

if (NOT APFEL_LIBRARIES)
  MESSAGE(STATUS "APFEL not found.")
  MESSAGE(STATUS "You can pass -DAPFEL_PREFIX=/path/to/apfel to cmake.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(APFEL DEFAULT_MSG APFEL_LIBRARIES APFEL_INCLUDE_DIRS)
