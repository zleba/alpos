pwd=$PWD
mkdir -p alposBuild
cd alposBuild
FNLODIR=$pwd/fastnloBuild
if [[ "$USER" = dnb ]]; then
    FNLODIR=/afs/ipp-garching.mpg.de/home/d/dnb/ATLAS/fastNLO/fastNLO/install
fi

cmakeCMD=cmake
if [ `command -v cmake3`  ]; then
    cmakeCMD=cmake3
fi
echo using $cmakeCMD

$cmakeCMD  ../alpos -DFNLO_PREFIX=$FNLODIR -DQCDNUM_PREFIX=$pwd/qcdnumBuild -DAPFEL_PREFIX=$pwd/apfel -DAPFELxx_PREFIX=$pwd/apfelxx -DCMAKE_INSTALL_PREFIX=$pwd/alposBuild

make -j`nproc`
#make install
