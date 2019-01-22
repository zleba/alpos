# - Try to find fastNLO
# Once done, this will define
#
# APPLGRID_FOUND
# APPLGRID_INCLUDE_DIRS
# APPLGRID_LIBRARIES

MESSAGE(STATUS "Find APPLgrid...")

FIND_LIBRARY(APPLGRID_LIBRARIES
  NAMES APPLgrid
  PATHS
  ${APPLGRID_PREFIX}/lib
  ${PREFIX}/lib
  $ENV{HOME}/.local/lib
  /usr/local/lib
  /usr/lib
  )

FIND_PATH(APPLGRID_INCLUDE_DIRS
  NAMES appl_grid/appl_grid.h
  PATHS
  ${APPLGRID_PREFIX}/include
  $ENV{HOME}/.local/include
  /usr/local/include
  /usr/include
  )

if (NOT APPLGRID_LIBRARIES)
  MESSAGE(STATUS "APPLGRID not found.")
  MESSAGE(STATUS "You can pass -DAPPLGRID_PREFIX=/path/to/applgrid to cmake.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(APPLgrid DEFAULT_MSG APPLGRID_LIBRARIES APPLGRID_INCLUDE_DIRS)
