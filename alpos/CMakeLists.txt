cmake_minimum_required (VERSION 2.6)
project (Alpos C Fortran CXX)
set (Alpos_VERSION_MAJOR 1)
set (Alpos_VERSION_MINOR 0)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

# use old ABI when compiling with gcc>=5
#add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)

# -- options for building Alpos

option(PYEXT
       "build Python bindings (requires CMake >= 3.8, Python and SWIG)"
       OFF)

# -- components distributed with Alpos

# AEM
file(GLOB aem_SRC "${PROJECT_SOURCE_DIR}/src/aem/*.f")
add_library(aem SHARED ${aem_SRC})

# APCALC
file(GLOB apc_SRC "${PROJECT_SOURCE_DIR}/src/apcalc/*.F")
add_library(apc SHARED ${apc_SRC})

# APCCPP
file(GLOB apccpp_SRC "${PROJECT_SOURCE_DIR}/src/apccpp/*.cc")
add_library(apccpp SHARED ${apccpp_SRC})
include_directories(SYSTEM ${PROJECT_SOURCE_DIR}/src/apccpp)

# specify component libraries for linking
set(FOUND_LIBRARIES aem apc apccpp gfortran)

# -- system components (required)
#find_package( Boost 1.66 COMPONENTS math REQUIRED )
#include_directories( ${Boost_INCLUDE_DIR} )



find_package(ROOT REQUIRED)
include_directories(SYSTEM ${ROOT_INCLUDE_DIRS})

find_package(LHAPDF REQUIRED)
include_directories(SYSTEM ${LHAPDF_INCLUDE_DIRS})

find_package(FastNLO REQUIRED)
include_directories(SYSTEM ${FNLO_INCLUDE_DIRS})

find_package(Eigen REQUIRED)
include_directories(SYSTEM ${EIGEN_INCLUDE_DIRS})

# append component libraries to variable used for linking
set(FOUND_LIBRARIES ${ROOT_LIBRARIES} ${LHAPDF_LIBRARIES} ${FNLO_LIBRARIES} ${FOUND_LIBRARIES})

# -- system components (optional)

# APFEL
find_package(APFEL)
if (APFEL_LIBRARIES)
    add_definitions(-D_CMAKE_FOUND_APFEL=1)
    include_directories(SYSTEM ${APFEL_INCLUDE_DIRS})
    set(FOUND_LIBRARIES ${FOUND_LIBRARIES} ${APFEL_LIBRARIES})
else()
    message(WARNING "Could not find APFEL: building Alpos without this interface")
endif()

# APFEL++
find_package(APFELxx)
if (APFELxx_LIBRARIES)
    add_definitions(-D_CMAKE_FOUND_APFELxx=1)
    include_directories(SYSTEM ${APFELxx_INCLUDE_DIRS})
    set(FOUND_LIBRARIES ${FOUND_LIBRARIES} ${APFELxx_LIBRARIES})
else()
    message(WARNING "Could not find APFEL++: building Alpos without this interface")
endif()

# QCDNUM
find_package(QCDNUM)
if (QCDNUM_LIBRARIES)
    add_definitions(-D_CMAKE_FOUND_QCDNUM=1)
    include_directories(SYSTEM ${QCDNUM_INCLUDE_DIRS})
    set(FOUND_LIBRARIES ${FOUND_LIBRARIES} ${QCDNUM_LIBRARIES})
else()
    message(WARNING "Could not find QCDNUM: building Alpos without this interface")
endif()

# APPLgrid
find_package(APPLgrid)
if (APPLGRID_LIBRARIES)
    add_definitions(-D_CMAKE_FOUND_APPLGRID=1)
    set(FOUND_LIBRARIES ${FOUND_LIBRARIES} ${APPLGRID_LIBRARIES_LIBRARIES})
else()
    message(WARNING "Could not find APPLgrid: building Alpos without this interface")
endif()

# -- target 'Alpos'

# choose which source files to compile
file(GLOB alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/*.cc" "${PROJECT_SOURCE_DIR}/src/alpos/tasks/*.cc" "${PROJECT_SOURCE_DIR}/src/alpos/functions/*.cc")

# specify 'include' directories
include_directories(BEFORE "${PROJECT_SOURCE_DIR}/include")

## exclude old source files
#list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AfastNLOnormDISalt.cc")
#list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AfastNLOalt.cc")

# exclude sources which require optional component 'QCDNUM'
if (NOT QCDNUM_LIBRARIES)
    #list(FILTER alpos_SRC EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/alpos/functions/AQcdnum(.*).cc")
    #list(FILTER alpos_SRC EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/alpos/functions/(.*)(n|N)ormDIS(alt)?.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AQcdnumAlphas.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AQcdnumDISCS.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AQcdnumDISCSEWFit.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AQcdnumInit.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AQcdnumPDF.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AfastNLOInterpolPDFasNormDIS.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AfastNLOnormDISalt.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AfastNLOnormDIS.cc")
endif()

# exclude sources which require optional component 'Apfel'
if (NOT APFEL_LIBRARIES)
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelAs.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelDISCS.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelDISCSEWFit.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelInit.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelPDF.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelQEDEvol.cc")
endif()

# exclude sources which require optional component 'Apfel++'
if (NOT APFELxx_LIBRARIES)
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelxxAlphas.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelxxDISCS.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApfelxxPDF.cc")
endif()

# exclude sources which require optional component 'APPLgrid'
if (NOT APPLGRID_LIBRARIES)
    #list(FILTER alpos_SRC EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApplgrid.cc")
    list(REMOVE_ITEM alpos_SRC "${PROJECT_SOURCE_DIR}/src/alpos/functions/AApplgrid.cc")
endif()

# build Alpos as a shared library
add_library(libalpos SHARED ${alpos_SRC})
set_target_properties(libalpos PROPERTIES OUTPUT_NAME alpos)

# specify the main executable
add_executable(alpos src/main.cc)

# need C++11 standard
set_property(TARGET alpos PROPERTY CXX_STANDARD 17)
set_property(TARGET alpos PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET libalpos PROPERTY CXX_STANDARD 17)
set_property(TARGET libalpos PROPERTY CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++17  -g")

# invoke linker from the C++ compiler
set_target_properties(alpos PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(libalpos PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(libalpos ${FOUND_LIBRARIES})
target_link_libraries(alpos libalpos)

if(PYEXT)
    cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
    # find SWIG package and load SWIG CMake macros
    find_package(SWIG REQUIRED)
    include(${SWIG_USE_FILE})

    # find Python libraries and include paths
    find_package(PythonLibs)
    include_directories(${PYTHON_INCLUDE_PATH})

    # specify SWIG source file for alpos
    set(alpos_SWIG_src ${PROJECT_SOURCE_DIR}/python/alpos.i)
    set_property(SOURCE ${alpos_SWIG_src} PROPERTY CPLUSPLUS ON)

    # add SWIG library (SWIG target name is '_alpos' by default)
    swig_add_library(alpos
                     LANGUAGE python
                     SOURCES ${alpos_SWIG_src})

    # remember to link SWIG wrapper to Alpos library
    target_link_libraries(_alpos libalpos)
endif(PYEXT)

# Install
install (TARGETS alpos
    DESTINATION bin)
install (TARGETS libalpos
    DESTINATION lib)
install (TARGETS aem
    DESTINATION lib)
install (TARGETS apc
    DESTINATION lib)
install (TARGETS apccpp
    DESTINATION lib)
install (DIRECTORY "${PROJECT_SOURCE_DIR}/include/alpos" DESTINATION include
    FILES_MATCHING PATTERN "*.h")
