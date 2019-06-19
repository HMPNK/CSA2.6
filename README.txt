# Chromosome Scale Assembler: 
# A high-throughput chromosome scale genome assembly pipeline for vertebrate genomes.

# Genome assembly of vertebrate genomes has improved much since long read technologies 
# approach read lengths larger than most repetitive elements.

# Yet, the final goal to achieve assembly of full chromosomes by sequencing data alone 
# has not been achieved.

# We have build a computationally highly efficient pipeline, which after contig assembly 
# performs post assembly improvements by ordering the assembly and closing gaps, as well 
# as splitting of low supported regions.

# The pipeline can use information from scaffolded assemblies (for example from HiC or
# 10X Genomics), or even from diverged (~65-100 Mya) reference genomes for ordering the
# contigs and thus support the assembly process. This typically results in improved 
# contig N50 when compared to current state of the art methods.

# For smaller vertebrate genomes (~1 Gbp) chromosome scale assemblies can be achieved 
# within 12h on high-end Desktop computers (Intel i7, 12 CPU threads, 128 GB RAM). 
# Larger mammalian genomes (~3Gbp) can be processed within 15-18 h on server equipment 
# (Xeon, 96 CPU threads, 1TB RAM).

# CSA can calculate (near) chromosomal scale assemblies from long read data and 
# comparisons with publicly available diverged reference genomes for Fish, Birds and Mammals.

![CSA Pipeline](https://github.com/HMPNK/CSA2.6/blob/master/Fig1.jpg)

# CSA-PIPELINE INSTALLATION
#
# THIRD PARTY OPEN SOURCE SOFTWARE USED (will be installed automatically by INSTALL.bash):
#
# Miniconda, wtdbg2.2, RAGOUT v1.0, Minimap2, last aligner v941, bedtools v2.27.1,
# samtools v1.7, pigz, GNU parallel, seqtk, seqkit, hgFakeAgp
#
# We deeply acknowledge the authors of these great tools! 

git clone https://github.com/HMPNK/CSA2.6.git

# INSTALL binaries and scripts
# !!! Make sure you have gcc, g++, make and zlib1g-dev installed on your system !!!
# Otherwise INSTALL.bash will fail in compilation steps
# On Ubuntu run:             "sudo apt install gcc g++ make zlib1g-dev"
# On Red Hat /CentOS 7 run:  "sudo yum install gcc gcc-c++ make zlib-devel"
# On OpenSuse run:           "sudo zypper install gcc gcc-c++ make zlib-devel"

cd CSA2.6
cd INSTALL
bash INSTALL.bash

# INITIALIZE MINICONDA ENVIRONMENT, ALTERNATIVELY YOU MAY JUST LOGOUT and LOGIN 
. ~/.bashrc

# RUN CSA2.6 PIPELINE

cd ../..
CSA2.6/CSA2.6c.pl

# INSTALLATION ISSUES

# We have tested CSA2.6c on fresh server installations of Red Hat 8 and Ubuntu 18.04/19.04,
# OpenSuse LEAP 15.1 and CentOS 7 as well as older Red Hat and Ubuntu versions.
#
# Our tests showed that RAGOUT V1.0 compilation needs gcc in version 4.6.2 or above.
# If you encounter the following error during the processing
#
# "perl: symbol lookup error: /tools/perl5/modules/lib/perl5/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined
# symbol: Perl_Istack_sp_ptr" 
#
# you have to add the correct perl lib location to your PERL5LIB environment variable
#
# export  PERL5LIB=/conda_root_directory_usually_homedir/.conda/pkgs/perl-5.22.0.1-0/lib/perl5:$PERL5LIB
#
# conda_root_directory_usually_homedir : has to be replaced with the correct location.

#TEST RUN ON YEAST OXFORD NANOPORE dataset:

mkdir CSA-TEST
cd CSA-TEST

#DOWNLOAD TEST DATA

wget https://nimbus.igb-berlin.de/index.php/s/njM2jqplusn17OZ/download
mv download SRR6476833.fa.gz
wget https://nimbus.igb-berlin.de/index.php/s/4VekUKms8tdL4V4/download
mv download sacCer.fa.gz

#CREATE CSA-PIPELINE SCRIPT

../CSA2.6/CSA2.6c.pl -r SRR6476833.fa.gz -g sacCer.fa.gz -t 4 -o SC_CSA -d SC_CSA > RUN-CSA-TEST.bash

#RUN CSA-PIPELINE

bash RUN-CSA-TEST.bash

#The test run will take about 30-60 min (you can speed it up by increasing "-t 4" on systems with higher CPU number)
