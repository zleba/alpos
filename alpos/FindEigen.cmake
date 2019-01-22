# - Try to find Eigen library
# Once done, this will define
#
# EIGEN_INCLUDE_DIRS

MESSAGE(STATUS "Try finding Eigen...")

FIND_PATH(EIGEN_INCLUDE_DIRS
  Eigen/Core
  /usr/include
  #/usr/include/eigen3
  /usr/local/include
  #/usr/local/include/eigen3
  ${EIGEN_PREFIX}/include
  )

if (NOT EIGEN_INCLUDE_DIRS)

  MESSAGE(STATUS "Eigen was not found on local system. Try to download it.")

  set(EIGEN_URL "http://bitbucket.org/eigen/eigen/get/3.2.8.tar.gz")
  set(EXTERNAL_FOLDER "external")
  set(EIGEN_DOWNLOAD_PATH "eigen.tar.gz")
  set(EIGEN_EXTRACTED_FILE "external/eigen")

  if (NOT EXISTS ${EXTERNAL_FOLDER})
    file(MAKE_DIRECTORY ${EXTERNAL_FOLDER})
  endif()

  # download tarball from web if not existing
  if (NOT EXISTS "${PROJECT_SOURCE_DIR}/${EXTERNAL_FOLDER}/${EIGEN_DOWNLOAD_PATH}")
    file(DOWNLOAD "${EIGEN_URL}" "${PROJECT_SOURCE_DIR}/${EXTERNAL_FOLDER}/${EIGEN_DOWNLOAD_PATH}")
  endif()

  # unpack the tarball
  if (NOT EXISTS ${PROJECT_SOURCE_DIR}/${EXTERNAL_FOLDER}/eigen)
    file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/${EXTERNAL_FOLDER}/eigen)
  endif()
  execute_process(
    COMMAND tar xzf ${EIGEN_DOWNLOAD_PATH} --strip-components=1 -C eigen
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/${EXTERNAL_FOLDER}
    )

  FIND_PATH(EIGEN_INCLUDE_DIRS
    Eigen/Core
    external/eigen/
  )
endif()

if (NOT EIGEN_INCLUDE_DIRS)
  MESSAGE(STATUS "Eigen library could neither be found locally nor be downloaded from the web.")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EIGEN DEFAULT_MSG EIGEN_INCLUDE_DIRS)
