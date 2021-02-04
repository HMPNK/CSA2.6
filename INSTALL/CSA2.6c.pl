#!/usr/bin/env perl

#CSA version 2.6 (Chromosome Scale Assembler, employing synteny with diverged reference genomes)

#THIS SCRIPT RUNS FOUR STEP CSA2.6 PIPELINE

#INPUT DATA: pacbio or oxford nanopore long reads (provided as fa.gz)
#	     suitable references (scaffolded assembly, 10X Genomics assembly, diverged reference genomes)

  use strict;
  use Getopt::Std;

my $userName =  $ENV{'LOGNAME'};

my $dir ="/home/$userName/CSA2.6";
my $wtdbg = "$dir/bin/wtdbg2.2";
my $script = "$dir/script";
my $bin = "$dir/bin";
my $reads = "";
my $out = "";
my $rundir = "";
my $ref = "";
my $threads = 2;
my $kmer = 21;
my $step = 4;
my $edge = 4;
my $ctg;
my $scaff;
my $guess = "on";
my $minmatch = 0.05;
my $cons = 0;
my $special = "-L 5000 -A";

#______________________________________________________________________
#get options:
my %options=();
getopts("C:S:o:r:g:t:k:s:e:d:x:m:p:l:", \%options);

if(! $ARGV[0] && ! $options{r})       {
       print STDERR "
CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)
Author: Heiner Kuhl, Phd (kuhl\@igb-berlin.de)

THIS SCRIPT RUNS FOUR STEP CSA2 PIPELINE

INPUT DATA: pacbio or oxford nanopore long reads (provided as fa.gz)
            suitable references (scaffolded assembly, 10X Genomics assembly, diverged reference genomes)

      Options:

              -d directory to run CSA in (will be created)

              -r input long reads as fa.gz (required!)

              -C input genome assembly contigs (overrides initial WTDBG assembly step
                                                fasta headers SHOULD NOT contain \":\", \"-\" or \".\" !
                                                Use simple fasta headers like \"scf1, scf2, ..., scfn\"
                                                Otherwise crashes STEP2 of CSA!)

              -g reference sequences (suitable backbone assembly or reference genomes,
                                     ordered by evolutionary distance (most similar first), 
                                     as comma separated list. Input format fasta or fasta.gz)

              -o output-prefix

              -t number threads (default=2, use multiples of 2)

              -k kmer wtdbg assembler (default=21, use 21 for high coverage data)

              -s step wtdbg assembler (default=4, use 4 for high coverage data and/or to save RAM)

              -e minimum edge coverage wtdbg assembler (default=4, decrease for low coverage data) 

              -m Min similarity in WTDBG (default=0.05)

              -x off/on (default on). Set \"-x on\" to further improve chromosomal assembly on the cost of more structural missassemblies.
                 Use \"-x on\" if you a have a closely related high quality reference genome with similar karyotype (family or genus level, same haploid chr no.)
                 If \"-x off\" a second assembly using guessed edges in the last assembly step will be output anyway since version 2.5 (\"xyz.final.B.fa\")

              -p 0,1,2 (0 = wtdbg-cons (default), 1 = wtpoa-cons; consensus calculation in step 1, 2 = wtdbg-cons -S 0
                        wtpoa-cns is slower but a bit more accurate; Parameter -S 0 is more stable on very long reads (e.g. ONT Ultra Long Reads))

              -l special parameters for wtdbg2 (default= -l \"-L 5000 -A\"; e.g. use -l \"-L 70000 --aln-min-length 25000 --keep-multiple-alignment-parts 1 -A\" for ultra long reads)

This script generates a bash script for running the pipeline! Write script to file and run by: nohup bash <script> & !
";       };
die "ERROR: NEED a read file!!!\n" if(! $options{r});
die "ERROR: unknown option $ARGV[0]" if $ARGV[0];

$reads = $options{r} if($options{r});
$out = $options{o} if($options{o});
$ref = $options{g} if($options{g});
$threads = $options{t} if($options{t});
$kmer = $options{k} if($options{k});
$step = $options{s} if($options{s});
$edge = $options{e} if($options{e});
$rundir = $options{d} if($options{d});
$ctg = $options{C} if($options{C});
$scaff = $options{S} if($options{S});
$guess = $options{x} if($options{x});
$minmatch = $options{m} if($options{m});
$cons = $options{p} if($options{p});
$special = $options{l} if($options{l});

if($threads<2){$threads=2};

die "\nERROR: Invalid (-r) long read file $reads ...Exiting.\n\n" if(! -e $reads);
die "\nERROR: Invalid (-r) does not match '.fa.gz' !!! Make sure to provide reads in gzipped fasta format! ...Exiting.\n\n" if(substr($reads, -5) ne "fa.gz");
die "\nERROR: Provide an output prefix! ...Exiting.\n\n" if($out eq "");
die "\nERROR: Specify a directory name for this RUN (will be created)! ...Exiting.\n\n" if($rundir eq "");
die "\nERROR: Paths are not allowed for RUN directory (no '\/')!!! ...Exiting.\n\n" if (index($rundir, "\/") > -1);
die "\nERROR: Provide path to reference assemblies as commaseparated list (ref1.fa,ref2.fa,..,refn.fa)! ...Exiting.\n\n" if($ref eq "");
if($ctg){die "\nERROR: Invalid (-C) contig file $ctg ...Exiting.\n\n" if(! -e $ctg);};
if($scaff){die "\nERROR: Invalid (-C) scaffolds file $scaff ...Exiting.\n\n" if(! -e $scaff);};

my @refcnt = split(/,/ , $ref);
foreach my $n (@refcnt) { die "$n does not exist ...Exiting." if(!-e $n);};

$ref =~ s/^/..\/..\//g;
$ref =~ s/,/,..\/..\//g;
$ref =~ s/,..\/..\/$//g;

$out =~ s/\./_/g;
$rundir =~ s/\./_/g;

#CSA2.6 assembly pipeline:
my $COMMAND="#!/usr/bash\nset -e\nset -o pipefail\n\n";

$COMMAND="$COMMAND\n#RUN CSA2.6 assembly pipeline:\n\n";

$COMMAND="$COMMAND

export PATH=$bin/RAGOUT_V1.0/lib:\$PATH
export PYTHONPATH=$bin/RAGOUT_V1.0/lib

mkdir $rundir
cd $rundir
";

if(! $ctg) {if(! $scaff){
$COMMAND="$COMMAND
mkdir 01_WTDBG
cd 01_WTDBG
$script/01_ASSEMBLY.pl -r ../../$reads -o $out -t $threads -k $kmer -s $step -e $edge -m $minmatch -p $cons -l \"$special\" > STEP_01.bash
time bash STEP_01.bash
cd ..
";
}}
elsif ( $ctg ) {
$COMMAND="$COMMAND
mkdir 01_WTDBG
cd 01_WTDBG
ln -s ../../$ctg $out.step1.fa
echo INPUT CONTIG STATS:
perl $bin/seq_n50.pl $out.step1.fa
cd ..
";

}

if(! $scaff) {
$COMMAND="$COMMAND

mkdir 02_ORDERCTG
cd 02_ORDERCTG
$script/02_ORDERCONTIGS.pl -c ../01_WTDBG/$out.step1.fa -g $ref -o $out -t $threads -x $guess > STEP_02.bash
time bash STEP_02.bash
cd ..
";
}

elsif( $scaff ) {
$COMMAND="$COMMAND

mkdir 02_ORDERCTG
cd 02_ORDERCTG
ln -s ../../$scaff $out.step2.fa
cd ..
";

}

$COMMAND="$COMMAND

mkdir 03_REASSEMBLE
cd 03_REASSEMBLE
$script/03_REASSEMBLE.pl -r ../../$reads -a ../02_ORDERCTG/$out.step2.fa -t $threads -o $out > STEP_03.bash
time bash STEP_03.bash
cd ..

mkdir 04_ORDERCTG
cd 04_ORDERCTG
$script/04_ORDERCONTIGS.pl -c ../03_REASSEMBLE/$out.step3.fa -g $ref -o $out -t $threads -x $guess > STEP_04.bash
time bash STEP_04.bash
cd ..

";
$COMMAND="$COMMAND\n#END CSA2.6 assembly pipeline\n\n";
print "$COMMAND\n";
