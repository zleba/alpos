pwd=$PWD
mkdir -p apfelxx
cd apfelxx
git clone --depth 1 https://github.com/vbertone/apfelxx.git
cd apfelxx
echo $PWD
#echo "set_property(TARGET dis PROPERTY POSITION_INDEPENDENT_CODE ON) " >> apfelxx/src/dis/CMakeLists.txt
#cat apfelxx/src/dis/CMakeLists.txt
cmake -DCMAKE_INSTALL_PREFIX=$pwd/apfelxx/ .
make -j`nproc`
make install
