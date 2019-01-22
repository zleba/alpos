# - Try to find Apfel
# Once done, this will define
#
# APFELxx_FOUND
# APFELxx_INCLUDE_DIRS
# APFELxx_LIBRARIES

MESSAGE(STATUS "Find Apfel++...")

FIND_LIBRARY(APFELxx_LIBRARIES
  NAMES apfelxx
  PATHS
  ${APFELxx_PREFIX}/lib
  ${PREFIX}/lib
  /usr/local/lib
  /usr/lib
  )

FIND_PATH(APFELxx_INCLUDE_DIRS
  apfel/dglap.h
  /usr/include
  /usr/local/include
  ${APFELxx_PREFIX}/include
  )

if (NOT APFELxx_LIBRARIES)
  MESSAGE(STATUS "APFEL++ not found.")
  MESSAGE(STATUS "You can pass -DAPFELxx_PREFIX=/path/to/apfel to cmake.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(APFELxx DEFAULT_MSG APFELxx_LIBRARIES APFELxx_INCLUDE_DIRS)
