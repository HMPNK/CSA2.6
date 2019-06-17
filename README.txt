#CSA-PIPELINE INSTALLATION

#unpack CSA

tar xvf CSA2.6_xyz.tar.gz

cd CSA2.6

#INSTALL binaries and scripts
#!!! Make sure you have gcc, g++, make and zlib1g-dev installed on your system !!!
#Otherwise INSTALL.bash will fail in compilation steps
#On ubuntu run "sudo apt install gcc g++ make zlib1g-dev"
#On red hat run "yum install gcc gcc-c++ make zlib-devel"
 

cd INSTALL
bash INSTALL.bash

#RUN CSA2.6 PIPELINE

cd ../..
CSA2.6/CSA2.6.pl
