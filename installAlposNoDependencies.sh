pwd=$PWD
mkdir -p alposBuild
cd alposBuild

cmakeCMD=cmake
if [ `command -v cmake3`  ]; then
    cmakeCMD=cmake3
fi
echo using $cmakeCMD

$cmakeCMD  ../alpos -DCMAKE_INSTALL_PREFIX=$pwd/alposBuild

make -j`nproc`
#make install
