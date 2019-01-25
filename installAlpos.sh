pwd=$PWD
mkdir -p alposBuild
cd alposBuild
FNLODIR=$pwd/fastnloBuild
if [[ "$USER" = dnb ]]; then
    FNLODIR=/afs/ipp-garching.mpg.de/home/d/dnb/ATLAS/fastNLO/fastNLO/install
fi
cmake3 ../alpos -DFNLO_PREFIX=$FNLODIR -DQCDNUM_PREFIX=$pwd/qcdnumBuild -DAPFELxx_PREFIX=$pwd/apfelxx -DCMAKE_INSTALL_PREFIX=$pwd/alposBuild

make -j`nproc`
#make install
