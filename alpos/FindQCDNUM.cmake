# - Try to find QCDNUM
# Once done, this will define
#
# QCDNUM_FOUND
# QCDNUM_LIBRARIES

MESSAGE(STATUS "Try to find QCDNUM")
set(QCDNUM_PREFIX $ENV{QCDNUM_PREFIX})

set(qcdnumlibs QCDNUM)
#set(qcdnumlibs hqstf mbutil qcdnum zmstf)
set(QCDNUM_LIBRARIES)

find_library(qcdnum_LIBRARY
  NAMES QCDNUM
  PATHS ${QCDNUM_PREFIX}/lib ${PREFIX}/lib /usr/lib /usr/local/lib
  )

FIND_PATH(QCDNUM_INCLUDE_DIRS
    NAMES QCDNUM/QCDNUM.h
  PATHS
  ${QCDNUM_PREFIX}/include
  $ENV{HOME}/.local/include
  /usr/local/include
  /usr/include
  )



list(APPEND QCDNUM_LIBRARIES ${qcdnum_LIBRARY})

#find_library(qcdnum_LIBRARY
#  NAMES qcdnum
#  PATHS ${QCDNUM_PREFIX}/lib ${PREFIX}/lib /usr/lib /usr/local/lib
#  )
#list(APPEND QCDNUM_LIBRARIES ${qcdnum_LIBRARY})

#find_library(hqstf_LIBRARY
#  NAMES hqstf
#  PATHS ${QCDNUM_PREFIX}/lib ${PREFIX}/lib /usr/lib /usr/local/lib
#  )
#list(APPEND QCDNUM_LIBRARIES ${hqstf_LIBRARY})

#find_library(mbutil_LIBRARY
#  NAMES mbutil
#  PATHS ${QCDNUM_PREFIX}/lib ${PREFIX}/lib /usr/lib /usr/local/lib
#  )
#list(APPEND QCDNUM_LIBRARIES ${mbutil_LIBRARY})

#find_library(zmstf_LIBRARY
#  NAMES zmstf
#  PATHS ${QCDNUM_PREFIX}/lib ${PREFIX}/lib /usr/lib /usr/local/lib
#  )
#list(APPEND QCDNUM_LIBRARIES ${zmstf_LIBRARY})

if (NOT QCDNUM_LIBRARIES)
  MESSAGE(STATUS "QCDNUM libraries not found.")
  MESSAGE(STATUS "You can pass -DQCDNUM_PREFIX=/path/to/qcdnum to cmake.")
endif()

find_package_handle_standard_args(QCDNUM DEFAULT_MSG QCDNUM_LIBRARIES)
