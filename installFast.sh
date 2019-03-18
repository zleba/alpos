#wget -qO-  https://fastnlo.hepforge.org/code/v23/fastnlo_toolkit-2.3.1pre-2411.tar.gz | tar zxv
#mkdir -p fastnloBuild
#mv fastnlo_toolkit-2.3.1pre-2411 fastnlo
rm -fr fastnlo
svn export https://github.com/zleba/diffPlotter/trunk/fastnlo
#exit

pwd=$PWD
cd fastnlo
./configure --prefix=$pwd/fastnloBuild  --disable-doxygen-doc  --with-lhapdf=$(dirname $(dirname `which lhapdf-config`)) CXXFLAGS=" -O3 -g -fopenmp "
make -j`nproc`
make install
