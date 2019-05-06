pwd=$PWD
mkdir -p alposBuild
cd alposBuild

CC=`which gcc`
CXX=`which g++`
FC=`which gfortran`
COMPILER_PATH=${CC/"\/bin\/gcc"/}

cmakeCMD=cmake
if [ `command -v cmake3`  ]; then
    cmakeCMD=cmake3
fi
echo using $cmakeCMD

$cmakeCMD  ../alpos -DCMAKE_INSTALL_PREFIX=$pwd/alposBuild # -DLHAPDF_PREFIX=`lhapdf-config --prefix`

make -j`nproc`
#make install
