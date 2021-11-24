cd pbsim2
./configure --prefix=`pwd`
make
make check
make install

cd ../minimap2 && make

cd ..
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
rm bedtools-2.29.1.tar.gz
cd bedtools2
make

cd ..
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar -xvf samtools-1.14.tar.bz2
rm samtools-1.14.tar.bz2
cd samtools-1.14   
./configure --prefix=`pwd`
make
make install
