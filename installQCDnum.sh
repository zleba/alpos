wget -qO-  https://www.nikhef.nl/~h24/qcdnum-files/download/qcdnum170114.tar.gz | tar zxv
mkdir -p qcdnumBuild
pwd=$PWD
mv qcdnum-17-01-14  qcdnum
cd qcdnum
./configure --prefix=$pwd/qcdnumBuild

sed -i 's/mxx0  = 300/mxx0  = 3000/;s/mqq0  = 150/mqq0  = 1500/;s/nwf0  = 1200000/nwf0  = 120000000/' qcdnum/inc/qcdnum.inc
sed -i 's/nhqstor = 6000000/nhqstor = 120000000/' hqstf/inc/hqstf.inc

make -j`nproc`
make install
