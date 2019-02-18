pwd=$PWD
mkdir -p apfel
cd apfel
git clone  --depth 1  https://github.com/scarrazza/apfel.git
cd apfel
echo $PWD

sed -i 's/parameter(nint_max=200)/parameter(nint_max=301)/;s/parameter(nint_max_DIS=120)/parameter(nint_max_DIS=220)/'   src/commons/grid.h

#echo "set_property(TARGET dis PROPERTY POSITION_INDEPENDENT_CODE ON) " >> apfelxx/src/dis/CMakeLists.txt
#cat apfelxx/src/dis/CMakeLists.txt
./configure --prefix=$pwd/apfel/ 
make -j`nproc`
make install
