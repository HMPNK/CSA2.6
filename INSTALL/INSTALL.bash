#!/usr/bash
set -e
set -o pipefail

##INSTALLS SOFTWARE USED IN CSA2.6 PIPELINE
rm -rf bin ../bin ../script ../CSA2.6.pl
mkdir bin

##INSTALL MINICONDA, if necessary

if hash conda 2>/dev/null; then
        echo "Conda is already available on your system. Proceeding..."
    else
        echo "conda is not available. First installing miniconda! Install to your home directory: $HOME ! Always answer \"yes\""
        sh ./Miniconda2-latest-Linux-x86_64.sh
        . ~/.bashrc
        echo "Conda is now available on your system. You might have to logout,login and restart INSTALL.bash. Trying to proceed..."
    fi


#UPDATE BASH TO ADD MINICONDA TO PATH
. ~/.bashrc

##MINIMAP2
conda install -y -c bioconda minimap2
TESTVAR=$(which minimap2)
ln -s $TESTVAR bin/

##LAST ALIGNER
conda install -y -c bioconda last=941
TESTVAR=$(which lastal)
ln -s $TESTVAR bin/
TESTVAR=$(which lastdb)
ln -s $TESTVAR bin/
TESTVAR=$(which last-split)
ln -s $TESTVAR bin/
TESTVAR=$(which maf-swap)
ln -s $TESTVAR bin/

##PARALLEL GZIP
conda install -y -c bioconda pigz
TESTVAR=$(which pigz)
ln -s $TESTVAR bin/

##SAMTOOLS
conda install -y -c bioconda samtools=1.7
TESTVAR=$(which samtools)
ln -s $TESTVAR bin/

##BEDTOOLS2
conda install -y -c bioconda bedtools=2.27.1
TESTVAR=$(which bedtools)
ln -s $TESTVAR bin/

##SEQTK
conda install -y -c bioconda seqtk
TESTVAR=$(which seqtk)
ln -s $TESTVAR bin/

##SEQKIT
conda install -y -c bioconda seqkit
TESTVAR=$(which seqkit)
ln -s $TESTVAR bin/

##HGFAKEAGP
conda install -y -c bioconda ucsc-hgfakeagp
TESTVAR=$(which hgFakeAgp)
ln -s $TESTVAR bin/

##GNU PARALLEL
conda install -y -c bioconda parallel
TESTVAR=$(which parallel)
ln -s $TESTVAR bin/

conda install -y -c bioconda openssl=1.0

##WTDBG2-2 (use sources as we need version 2.2!)
tar xvf wtdbg-v2.2.tar.gz
cd wtdbg2-2.2
make
rm *.h *.c
cd ..
mv wtdbg2-2.2 bin/wtdbg2.2

##RAGOUT_V1.0 (old version performs better, if it does not compile try "module load compiler/gcc-7.2.0")
tar xvf RAGOUT_V1.0.tar.gz
cd RAGOUT_V1.0
make clean
make
cd ..
mv RAGOUT_V1.0 bin/

##PUT BINARIES AND SCRIPTS IN THE RIGHT PLACE
mv bin ../bin
cp -r ./script ../script

##CHANGE CSA2.6 directory path in scripts
REPLACE=/home/\$userName/CSA2.6
NEW=$PWD/..
sed "s|$REPLACE|$NEW|g" ./CSA2.6c.pl > ../CSA2.6c.pl
sed "s|$REPLACE|$NEW|g" ./01_ASSEMBLY.pl > ../script/01_ASSEMBLY.pl
sed "s|$REPLACE|$NEW|g" ./02_ORDERCONTIGS.pl > ../script/02_ORDERCONTIGS.pl
sed "s|$REPLACE|$NEW|g" ./03_REASSEMBLE.pl > ../script/03_REASSEMBLE.pl
sed "s|$REPLACE|$NEW|g" ./04_ORDERCONTIGS.pl > ../script/04_ORDERCONTIGS.pl

chmod 770 ../CSA2.6c.pl ../script/01_ASSEMBLY.pl ../script/02_ORDERCONTIGS.pl ../script/03_REASSEMBLE.pl ../script/04_ORDERCONTIGS.pl

#TEST binaries
#TEST binaries
../bin/bedtools -h
printf ">1\nACGT" | ../bin/hgFakeAgp stdin stdout
../bin/lastal -h
../bin/lastdb -h
../bin/last-split -h
../bin/maf-swap -h
../bin/minimap2 -h
printf "echo ABCD" | ../bin/parallel
../bin/pigz -h
../bin/samtools sort
../bin/seqkit -h
echo "ACGT" | ../bin/seqtk seq -
../bin/wtdbg2.2/kbm2 -V
../bin/wtdbg2.2/pgzf -V
../bin/wtdbg2.2/wtdbg2 -V
../bin/wtdbg2.2/wtdbg-cns -V
../bin/wtdbg2.2/wtpoa-cns -V
../bin/RAGOUT_V1.0/ragout.py -h
touch test.maf
../bin/RAGOUT_V1.0/lib/ragout-maf2synteny test.maf MAF test.maf 100 1000
rm test.maf MAF -rf

##FINISH
echo "#######################################"
echo "# CSA2.6 INSTALLATION was successful! #"
echo "#######################################"
