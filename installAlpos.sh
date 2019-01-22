pwd=$PWD
mkdir -p alposBuild
cd alposBuild
cmake ../alpos -DFNLO_PREFIX=$pwd/fastnloBuild -DQCDNUM_PREFIX=$pwd/qcdnumBuild -DCMAKE_INSTALL_PREFIX=$pwd/alposBuild
make -j`nproc`
make install
