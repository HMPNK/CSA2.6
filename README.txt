# CSA-PIPELINE INSTALLATION

# 

git clone https://github.com/HMPNK/CSA2.6.git

cd CSA2.6

# INSTALL binaries and scripts
# !!! Make sure you have gcc, g++, make and zlib1g-dev installed on your system !!!
# Otherwise INSTALL.bash will fail in compilation steps
# On ubuntu run "sudo apt install gcc g++ make zlib1g-dev"
# On red hat run "yum install gcc gcc-c++ make zlib-devel"
 

cd INSTALL
bash INSTALL.bash

# RUN CSA2.6 PIPELINE

cd ../..
CSA2.6/CSA2.6c.pl


# INSTALLATION ISSUES

# We have tested CSA2.6c on fresh server installations of Red Hat 8 and Ubuntu 18.04/19.04 as well as older Red Hat and Ubuntu versions.
#
# Our tests showed that RAGOUT V1.0 compilation needs gcc in version 4.6.2 or above.
# If you encounter the following error during the processing
#
# "perl: symbol lookup error: /tools/perl5/modules/lib/perl5/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined symbol: Perl_Istack_sp_ptr" 
#
# you have to add the correct perl lib location to your PERL5LIB environment variable
#
# export  PERL5LIB=/conde_root_directory_usually_homedir/.conda/pkgs/perl-5.22.0.1-0/lib/perl5:$PERL5LIB
#
# conde_root_directory_usually_homedir : has to be replace with the correct location.
