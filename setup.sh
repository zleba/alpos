if [ -d "/cvmfs/sft.cern.ch/lcg" ]; then
source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-slc6-gcc7-opt/setup.sh
export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current/:`lhapdf-config --datadir`
else
echo "The cvmfs file system not mounted"
fi

export PATH=$PWD/qcdnumBuild/bin:$PATH
export PATH=$PWD/fastnloBuild/bin:$PATH
export PATH=$PWD/alposBuild:$PATH
export ALPOS_DIR=$PWD/alpos
