#!/usr/bash
set -e
set -o pipefail

if [ -f RAGOUT_V1.0.tar.gz ]; then
    echo "INSTALLING CSA2.6c"
    else
    echo "You have to be in the directory where INSTALL.bash is located!!!"
    exit;
fi

##INSTALLS SOFTWARE USED IN CSA2.6 PIPELINE
rm -rf bin ../bin ../script ../CSA2.6.pl
mkdir bin

##INSTALL MINICONDA, if necessary

if hash conda 2>/dev/null; then
        echo "Conda is already available on your system. Proceeding..."
    else
        echo "conda is not available. First installing miniconda! Install to your home directory: $HOME ! Always answer \"yes\""
        wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh
        sh ./Miniconda2-latest-Linux-x86_64.sh
        source ~/.bashrc
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

#conda install -y -c bioconda openssl=1.0 #version should not be fixed anymore 24/10/10
conda install -y -c bioconda openssl

##WTDBG2-2 (use sources as we need version 2.2!)
tar xvf wtdbg-v2.2.tar.gz
cd wtdbg2-2.2
make
rm *.h *.c
cd ..
mv wtdbg2-2.2 bin/wtdbg2.2
ln -s ./wtdbg2.2/scripts/seq_n50.pl bin/seq_n50.pl

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
echo "Testing binaries..........."
echo ""

if ../bin/bedtools -h >/dev/null; then
        echo "bedtools.................ok"
    else
        echo "bedtools.............failed"
fi

if printf ">1\nACGT" | ../bin/hgFakeAgp stdin stdout >/dev/null; then
        echo "hgFakeAgp................ok"
    else
        echo "hgFakeAgp............failed"
fi

if ../bin/lastal -h >/dev/null; then
        echo "lastal...................ok"
    else
        echo "lastal...............failed"
fi

if ../bin/lastdb -h >/dev/null; then
        echo "lastdb...................ok"
    else
        echo "lastdb...............failed"
fi

if ../bin/last-split -h >/dev/null; then
        echo "last-split...............ok"
    else
        echo "last-split...........failed"
fi

if ../bin/maf-swap -h >/dev/null; then
        echo "maf-swap.................ok"
    else
        echo "maf-swap.............failed"
fi

if ../bin/minimap2 -h >/dev/null; then
        echo "minimap2.................ok"
    else
        echo "minimap2.............failed"
fi

if printf "echo ABCD" | ../bin/parallel >/dev/null 2>&1; then
        echo "parallel.................ok"
    else
        echo "parallel.............failed"
fi

if ../bin/pigz -h >/dev/null 2>&1; then
        echo "pigz.....................ok"
    else
        echo "pigz.................failed"
fi

if printf ">1\nacgt" | ../bin/samtools faidx -;rm ./-.fai; then
        echo "samtools.................ok"
    else
        echo "samtools.............failed"
fi

if ../bin/seqkit -h >/dev/null; then
        echo "seqkit...................ok"
    else
        echo "seqkit...............failed"
fi

if printf ">1\nACGT" | ../bin/seqtk seq - >/dev/null; then
        echo "seqtk....................ok"
    else
        echo "seqtk................failed"
fi

if ../bin/wtdbg2.2/kbm2 -V >/dev/null; then
        echo "kbm2.....................ok"
    else
        echo "kbm2.................failed"
fi

if ../bin/wtdbg2.2/pgzf -V >/dev/null; then
        echo "pgzf.....................ok"
    else
        echo "pgzf.................failed"
fi

if ../bin/wtdbg2.2/wtdbg2 -V >/dev/null; then
        echo "wtdbg2...................ok"
    else
        echo "wtdbg2...............failed"
fi

if ../bin/wtdbg2.2/wtdbg-cns -V >/dev/null; then
        echo "wtdbg-cns................ok"
    else
        echo "wtdbg-cns............failed"
fi

if ../bin/wtdbg2.2/wtpoa-cns -V >/dev/null; then
        echo "wtpoa-cns................ok"
    else
        echo "wtpoa-cns............failed"
fi

if ../bin/RAGOUT_V1.0/ragout.py -h >/dev/null; then
        echo "ragout.py................ok"
    else
        echo "ragout.py............failed"
fi

touch test.maf
if ../bin/RAGOUT_V1.0/lib/ragout-maf2synteny test.maf MAF test.maf 100 1000 >/dev/null 2>&1; then
        echo "ragout-maf2synteny.......ok"
    else
        echo "ragout-maf2synteny...failed"
fi
rm test.maf MAF -rf

echo ""

##FINISH
echo "#######################################"
echo "# CSA2.6 INSTALLATION was successful! #"
echo "#######################################"
