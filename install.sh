cd pbsim2
./configure --prefix=`pwd`
make
make check
make install

cd ../minimap2 && make
