mkdir -p f2ave
cd f2ave
F2DIR=$PWD
echo $F2DIR
wget "https://svnsrv.desy.de/viewvc/heraverager/trunk/?view=tar" 
ls
mv "index.html?view=tar" heraverager-trunk.tar.gz
tar xzvf heraverager-trunk.tar.gz
cd trunk
autoreconf -i
./configure --prefix=$F2DIR
make -j`nproc`
make install
