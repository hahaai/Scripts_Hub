################ install singularity 
# https://singularity.lbl.gov/install-linux
# https://singularity.lbl.gov/docs-installation
sudo yum groupinstall "Development Tools"
sudo yum install libarchive-devel
sudo yum install squashfs-tools

VERSION=2.5.2
wget https://github.com/singularityware/singularity/releases/download/$VERSION/singularity-$VERSION.tar.gz
tar xvf singularity-$VERSION.tar.gz
cd singularity-$VERSION
./configure --prefix=/usr/local
make
sudo make install



################ parallel

wget https://ftpmirror.gnu.org/parallel/parallel-20200322.tar.bz2
wget https://ftpmirror.gnu.org/parallel/parallel-20200322.tar.bz2.sig
gpg parallel-20200322.tar.bz2.sig
bzip2 -dc parallel-20200322.tar.bz2 | tar xvf -
cd parallel-20200322
./configure && make && sudo make install