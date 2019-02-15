pwd=$PWD
mkdir -p apfel
cd apfel
git clone  https://github.com/scarrazza/apfel.git
cd apfel
echo $PWD
#echo "set_property(TARGET dis PROPERTY POSITION_INDEPENDENT_CODE ON) " >> apfelxx/src/dis/CMakeLists.txt
#cat apfelxx/src/dis/CMakeLists.txt
./configure --prefix=$pwd/apfel/ 
make -j`nproc`
make install
