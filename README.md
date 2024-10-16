# Chromosome Scale Assembler 
## A high-throughput chromosome scale genome assembly pipeline for vertebrate genomes.

Genome assembly of vertebrate genomes has improved much since long read technologies 
approach read lengths larger than most repetitive elements.

Yet, the final goal to achieve assembly of full chromosomes by sequencing data alone 
has not been achieved.

We have build a computationally highly efficient pipeline, which after contig assembly 
performs post assembly improvements by ordering the assembly and closing gaps, as well 
as splitting of low supported regions.

![CSA-PIPELINE](https://github.com/HMPNK/CSA2.6/blob/master/Fig1.png)

The pipeline can use information from scaffolded assemblies (for example from HiC or
10X Genomics), or even from diverged (~65-100 Mya) reference genomes for ordering the
contigs and thus support the assembly process. This typically results in improved 
contig N50 when compared to current state of the art methods.

For smaller vertebrate genomes (~1 Gbp) chromosome scale assemblies can be achieved 
within 12h on high-end Desktop computers (Intel i7, 12 CPU threads, 128 GB RAM). 
Larger mammalian genomes (~3Gbp) can be processed within 15-18 h on server equipment 
(Xeon, 96 CPU threads, 1TB RAM).

CSA can calculate (near) chromosomal scale assemblies from long read data and 
comparisons with publicly available diverged reference genomes for Fish, Birds and Mammals.

# CSA-PIPELINE INSTALLATION

 THIRD PARTY OPEN SOURCE SOFTWARE USED (will be installed automatically by INSTALL.bash):

 Miniconda, wtdbg2.2, RAGOUT v1.0, Minimap2, last aligner v941, bedtools v2.27.1,
 samtools v1.7, pigz, GNU parallel, seqtk, seqkit, hgFakeAgp

 We deeply acknowledge the authors of these great tools! 

```sh
git clone https://github.com/HMPNK/CSA2.6.git
```
 Make sure you have gcc, g++, make and zlib1g-dev installed on your system !!!
 Otherwise INSTALL.bash will fail in compilation steps.

 On Ubuntu run:
 ```sh
 sudo apt install gcc g++ make zlib1g-dev
 ```
  On Red Hat /CentOS 7 run:
```sh
sudo yum install gcc gcc-c++ make zlib-devel
```
 On OpenSuse run:
 ```sh
 sudo zypper install gcc gcc-c++ make zlib-devel
```
INSTALL CSA binaries and scripts!!! 
```sh
cd CSA2.6
cd INSTALL
bash INSTALL.bash
#alterantively try "INSTALL-Oct24.bash", if INSTALL.bash does not finish correctly.
```
 INITIALIZE MINICONDA ENVIRONMENT, ALTERNATIVELY YOU MAY JUST LOGOUT and LOGIN 
```sh
. ~/.bashrc
```
 RUN CSA2.6 PIPELINE
```sh
cd ../..
CSA2.6/CSA2.6c.pl
```
 INSTALLATION ISSUES

 We have tested CSA2.6c on fresh server installations of Red Hat 8 and Ubuntu 18.04/19.04,
 OpenSuse LEAP 15.1 and CentOS 7 as well as older Red Hat and Ubuntu versions.

 Our tests showed that RAGOUT V1.0 compilation needs gcc in version 4.6.2 or above.
 If you encounter the following error during the processing

 "perl: symbol lookup error: /tools/perl5/modules/lib/perl5/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined
 symbol: Perl_Istack_sp_ptr" 

 you have to add the correct perl lib location to your PERL5LIB environment variable
```sh
 export  PERL5LIB=/conda_root_directory_usually_homedir/.conda/pkgs/perl-5.22.0.1-0/lib/perl5:$PERL5LIB
```
 conda_root_directory_usually_homedir : has to be replaced with the correct location.

Some users encountered problems in Step3 of the pipeline, if the language variable was NOT set to "en_US.UTF-8".
Make sure that your system uses "en_US.UTF-8" by:
```sh
export LANG='en_US.UTF-8'
```
If you encounter problems during the run regarding paths to the longread-file or reference file locations. Try to link these files to the folder where you run CSA and omit paths in the CSA commandline.


# Default parameters
CSA default parameters are currently tweaked for Pacbio RSII and ONT reads (30-60X, N50 readlength 10-30 kbp)
We have found that some SEQUEL datasets behave quite different, here adding custom parameters 
for WTDBG2 will help: -l "-p 0 -k 15 -L5000 -S 2 -A" .

Ultra long read assembly may be improved by increasing read length cut-off and alignment length cut-off:
For example: -l "-L 70000 --aln-min-length 25000 --keep-multiple-alignment-parts 1 -A" , worked well for a human dataset with 
readlength N50 of 70 kbp.

# TEST RUN ON YEAST OXFORD NANOPORE dataset:
```sh
mkdir CSA-TEST
cd CSA-TEST
```
DOWNLOAD TEST DATA
```sh
wget https://nimbus.igb-berlin.de/index.php/s/njM2jqplusn17OZ/download
mv download SRR6476833.fa.gz
wget https://nimbus.igb-berlin.de/index.php/s/4VekUKms8tdL4V4/download
mv download sacCer.fa.gz
```
CREATE CSA-PIPELINE SCRIPT
```sh
../CSA2.6/CSA2.6c.pl -r SRR6476833.fa.gz -g sacCer.fa.gz -t 4 -o SC_CSA -d SC_CSA > RUN-CSA-TEST.bash
```
RUN CSA-PIPELINE
```sh
bash RUN-CSA-TEST.bash
```
The test run will take about 30-60 min (you can speed it up by increasing "-t 4" on systems with higher CPU number)

You may see a few "core dumps", these are single gap re-assembly jobs that failed, you should just not care. Some gap re-assembly jobs may run very long due to huge read pile ups in that location, therefore we force gap re-assemblies to stop after 20 minutes using "timeout".

# POLISHING CONSENSUS QUALITY OF CSA ASSEMBLIES:
All of the widely used genome assembly tools like CANU, FALCON, FLYE and even SHASTA do not include the polishing algorithms that are needed to get error rates above Q40 (which is needed for annotation without an excess of frameshift errors). As it is currently common practice, we leave the choice of consensus polishing to the user. This is also, because choosing the right polishing pipeline for the different long-read technologies (SMRT or ONT) can be a complex topic (e. g. involving different sequencing libraries and also depending on properties of the species you are sequencing) and in our opinion requires some human decision making.

We have currently made good experience doing one iteration of polishing by the FLYE polisher and long read data, followed by two iterations using PILON and illumina short read data. Some long read polishers (e.g. MEDAKA and FLYE) do not like "n" characters, it is necessary to split the CSA scaffolds into contigs before polishing and re-build the scaffolds afterwards.


# FUNDING
This work was funded by the German Research foundation (DFG) “eigene Stelle” grant within the project KU 3596/1-1; project number: 324050651.

# CITE
Heiner Kuhl, Ling Li, Sven Wuertz, Matthias Stöck, Xu-Fang Liang, Christophe Klopp, CSA: A high-throughput chromosome-scale assembly pipeline for vertebrate genomes, GigaScience, Volume 9, Issue 5, May 2020, giaa034, https://doi.org/10.1093/gigascience/giaa034 
