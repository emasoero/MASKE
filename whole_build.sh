set -e
echo "Whole Build Process"

rootDir=$PWD

echo "Going to lammps..."

cd $rootDir/lammps
echo "Removing lammps/build folder..."
rm -rf build
echo "Making lammps/build folder..."
mkdir build
cd build
echo $PWD

echo "Running cmake command for lammps..."

cmake -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_BUILD_TYPE=Release -DWITH_JPEG=Off -DWITH_PNG=Off -DWITH_FFMPEG=Off -DBUILD_LIB=On -DBUILD_OMP=Off -DPKG_MISC=yes -DPKG_USER-MASKE=yes ../cmake

echo "Running make command"
make -j2

echo "lammps building completed..."

echo "Going to root directory..."
cd $rootDir
echo $PWD

echo "Remove build folder..."
rm -rf build

echo "Making build folder..."
mkdir build
cd build

echo $PWD

set +e

echo "Running cmake command for build folder..."
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_NUFEB=On ..

echo "Going to NUFEB folder..."
cd $rootDir/NUFEB
echo $PWD

./uninstall.sh
echo "Uninstalling existing NUFEB"

echo "Running install.sh script in NUFEB"
./install.sh --enable-misc --static

set -e

echo "Going to build folder again..."
cd $rootDir/build

echo $PWD

echo "Running cmake command again..."
cmake -DCMAKE_CXX_FLAGS="-std=c++11" -DCMAKE_BUILD_TYPE=Release -DWITH_NUFEB=On ..

echo "Running final make command..."
make -j2
