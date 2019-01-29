wget -qO-  https://fastnlo.hepforge.org/code/v23/fastnlo_toolkit-2.3.1pre-2411.tar.gz | tar zxv
mkdir -p fastnloBuild
mv fastnlo_toolkit-2.3.1pre-2411 fastnlo
pwd=$PWD
cd fastnlo
./configure --prefix=$pwd/fastnloBuild  --disable-doxygen-doc  --with-lhapdf=$(dirname $(dirname `which lhapdf-config`))
make -j`nproc`
make install
