#!/usr/bin/env perl

#CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)

#THIS SCRIPT RUNS WTDBG assembly as first step of the pipeline

#INPUT DATA: pacbio or oxford nanopore long reads (provided as fa.gz)


  use strict;
  use Getopt::Std;
my $userName =  $ENV{'LOGNAME'};

my $dir ="/home/$userName/CSA2.6";
my $wtdbg = "$dir/bin/wtdbg2.2";
my $bin = "$dir/bin";
my $script = "$dir/script";
my $reads = "";
my $out = "CSA";
my $ref = "";
my $threads = 1;
my $kmer = 21;
my $step = 4;
my $edge = 4;
my $minmatch = 0.05;
my $cons = "wtdbg-cns";
my $special = "-L 5000 -A";
my $lrpe = 50000;

#______________________________________________________________________
#get options:
my %options=();
getopts("o:r:g:t:k:s:e:m:p:l:", \%options);

if(! $ARGV[0] && ! $options{r})       {
       print STDERR "
CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)
Author: Heiner Kuhl, Phd (kuhl\@igb-berlin.de)

THIS SCRIPT RUNS WTDBG assembly as first step of the pipeline

INPUT DATA: pacbio or oxford nanopore long reads (provided as fa.gz)


      Options:\n
              -r input long reads as fa.gz (required!)
              -o output-prefix
              -t number threads (default=1)
              -k kmer wtdbg assembler (default=21)
              -s step wtdbg assembler (default=4)
              -e minimum edge coverage wtdbg assembler (default=4) 
	      -m Min similarity in WTDBG (default=0.05)
              -p 1 = wtpoa-cns 0 = wtdbg-cns 2= wtdbg-cns -S 0(default = 0) 
		 (wtpoa-cns longer runtime, slightly better consensus accuracy; Parameter -S 0 is more stable on ultra long reads)
              -l special parameters for wtdbg2 (default= -l \"-L 5000 -A\"; e.g. use -l \"-L 70000 --aln-min-length 25000 --keep-multiple-alignment-parts 1 -A\" for ultra long reads)
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
$minmatch = $options{m} if($options{m});
$cons = "wtpoa-cns" if($options{p}==1);
$cons = "wtdbg-cns -S 0" if($options{p}==2);
$special = $options{l} if($options{l});

if($options{p}==2) {$lrpe=150000};

my $threads2 =int($threads/2);

die "\nERROR: Invalid (-r) long read file $reads ...Exiting.\n\n" if(! -e $reads);
die "\nERROR: Invalid (-r) does not match '.fa.gz' !!! Make sure to provide reads in gzipped fasta format! ...Exiting.\n\n" if(substr($reads, -5) ne "fa.gz");

my $minmatch2 = $minmatch+0.025;

#01_WTDBG assembly:
my $COMMAND="#!/usr/bash\nset -e\nset -o pipefail\n\n";

$COMMAND="$COMMAND\n\necho;date;echo RUN WTDBG on $reads;echo\n\n";
$COMMAND="$COMMAND

#RUN TWO WTDBG ASSEMBLIES
$wtdbg/wtdbg2 -t $threads -i $reads -o $out-WTDBG -p $kmer -S $step -e $edge -s $minmatch $special 2>$out-wtdbg.log
$wtdbg/wtdbg2 -t $threads -i $reads -o $out-WTDBG_2 --load-alignments $out-WTDBG.alignments.gz -p $kmer -S $step -e $edge -s $minmatch2 $special 2>>$out-wtdbg.log

echo;date;echo CREATE CONSENSUS / use $cons;echo


printf \"$wtdbg/$cons -i $out-WTDBG.ctg.lay.gz -f -t $threads2 -o $out-wtdbg.fa 2>>$out-wtdbg.log \\n $wtdbg/$cons -i $out-WTDBG_2.ctg.lay.gz -f -t $threads2 -o $out-wtdbg_2.fa 2>>$out-wtdbg.log\" | $bin/parallel -j 2 > parallel.log 2>&1


echo;date;echo COMPARE ASSEMBLIES AND SPLIT CONTIGS AT DIFFERENCES;echo

$bin/minimap2 -I 100G -t $threads -x asm5 $out-wtdbg.fa $out-wtdbg_2.fa > $out-wtdbg_1vs2.paf 2>1vs2.minimap.log
awk '{if(\$12>=60 && \$9-\$8 > 30000){print \$1\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9}}' $out-wtdbg_1vs2.paf| sort -k2,2V -k3,3n -k4,4n| awk '{if(\$2==oq && \$1!=oR){n=split(ol,d,\"\\t\");if(d[n]<=\$3){print \$2\"\\t\"d[n]\"\\t\"\$3} else {print \$2\"\\t\"\$3\"\\t\"d[n]}};oR=\$1;oq=\$2;ol=\$0}' > $out-wtdbg.fa.SPLIT.bed

if [ -s $out-wtdbg.fa.SPLIT.bed ]
then
$bin/samtools faidx $out-wtdbg.fa
cut -f 1,2 $out-wtdbg.fa.fai > $out-wtdbg.fa.sizes
$bin/bedtools complement -g $out-wtdbg.fa.sizes -i $out-wtdbg.fa.SPLIT.bed > $out-wtdbg.fa.GOOD.bed
else
$bin/samtools faidx $out-wtdbg.fa
cut -f 1,2 $out-wtdbg.fa.fai > $out-wtdbg.fa.sizes
awk '{print \$1\"\\t0\\t\"\$2}' $out-wtdbg.fa.sizes > $out-wtdbg.fa.GOOD.bed
fi

sort -k1,1V -k2,2n $out-wtdbg.fa.SPLIT.bed $out-wtdbg.fa.GOOD.bed | $bin/bedtools merge -d -1000 | awk '{print \$0\"\\t\"\$1\"_\"x[\$1]++}' > $out-wtdbg.fa.SPLIT.REGIONS.bed
$bin/bedtools getfasta -name -fi $out-wtdbg.fa -fo $out.step1.fa -bed $out-wtdbg.fa.SPLIT.REGIONS.bed

##WORKARROUND TO ASSEMBLE MISSING REGIONS IN WTDBG2 ASSEMBLIES (TYPICAL 10-20 MB in a 1GB genome)
$bin/pigz -dc $reads | awk '{if(substr(\$1,1,1)==\">\") {print substr(\$1,2) > \"readnames.list\" ;print} else{print}}' | $bin/minimap2 -I 100G -t $threads -x map-pb $out.step1.fa - > $out.step1.fa.paf 2> remap.log
awk '{if(o!=\$1 && (\$4-\$3)/\$2>0.1) {print \$1};o=\$1}' $out.step1.fa.paf | sort -u > reads_mapped.list
sort reads_mapped.list readnames.list|uniq -u > reads_UNmapped.list
$bin/pigz -dc $reads | $bin/seqtk subseq /dev/stdin reads_UNmapped.list | $bin/pigz -c > reads_UNmapped.fa.gz
$wtdbg/wtdbg2 -t $threads -i reads_UNmapped.fa.gz -o UNmapped-wtdbg2 -p $kmer -S $step -e $edge -s $minmatch -L 1000 -A  2>UNmapped-wtdbg.log
$wtdbg/$cons -i UNmapped-wtdbg2.ctg.lay.gz -f -t 90 -o UNmapped-wtdbg2.fa 2>>UNmapped-wtdbg.log
$bin/seqtk seq -l 0 UNmapped-wtdbg2.fa >> $out.step1.fa

echo;date;echo ASSEMBLY IMPROVEMENTS BY READ RE-MAPPING AND SPLITTING;echo

#MAP LRs and LRPEs

printf \"$bin/minimap2 -I 100G -t $threads2 -x map-pb $out.step1.fa $reads > $out.step1.LR.paf 2>LR.log \\n $bin/pigz -dc $reads | $bin/seqtk seq -l 0 - | awk -f $script/EOF_LR.awk | $bin/minimap2 -I 100G -t $threads2 -x sr -k19 -w10 -a $out.step1.fa - 2>LRPE.log | $bin/samtools view -Sb -@ 6 - > $out.step1.LRPE.bam\" | $bin/parallel -j 2 >> parallel.log 2>&1


#Get genome coverage by best matches of LRs:
$bin/seqtk comp $out.step1.fa | cut -f 1,2 | sort -k1,1V -k2,2n > $out.step1.sizes
awk '{if(\$12>=60 && o!=\$1){print};o=\$1}' $out.step1.LR.paf | cut -f 6,8,9 | sort -k1,1V -k2,2n | $bin/bedtools genomecov -bga -i - -g $out.step1.sizes > $out.step1.LR.bedcov

#Get Genome coverage by good quality (MQ20) and concordantly mapped LRPEs:
$bin/samtools view -h $out.step1.LRPE.bam | awk '{if(substr(\$1,1,1)==\"@\" || \$2<1000){print}}'| $bin/samtools view -Sbu - | $bin/bedtools bamtobed -bedpe -i - | awk '{if(\$1==\$4 && \$8>=20 && \$9==\$10 && \$6-\$2<$lrpe){print}}'| cut -f 1,2,6| sort -k1,1V -k2,2n|$bin/bedtools genomecov -bga -i - -g $out.step1.sizes > $out.step1.LRPE.bedcov

#JOIN COVERAGE FILES
$bin/bedtools unionbedg -i $out.step1.LR.bedcov  $out.step1.LRPE.bedcov -g $out.step1.sizes | sort -k1,1V -k2,2n > $out.step1.LR_LRPE.bedcov

#Get good regions
awk '{if(\$4<1 && \$5<1) {print}}' $out.step1.LR_LRPE.bedcov| $bin/bedtools merge | $bin/bedtools complement -g $out.step1.sizes -i - | awk '{if(\$2<\$3){print}}' > $out.step1.LR_LRPE.SPLIT.bed
$bin/bedtools complement -g $out.step1.sizes -i $out.step1.LR_LRPE.SPLIT.bed | awk '{if(\$2<\$3){print}}' > $out.step1.LR_LRPE.CUT.bed
sort -k1,1V -k2,2n $out.step1.LR_LRPE.SPLIT.bed $out.step1.LR_LRPE.CUT.bed | awk  '{if(\$3>\$2 && \$3-\$2>2000){print \$0\"\\t\"\$1\"_\"d[\$1]++}}' > $out.step1.LR_LRPE.GOOD.bed

$bin/bedtools getfasta -bed $out.step1.LR_LRPE.GOOD.bed -fi $out.step1.fa -fo $out.step1b.fa -name > /dev/null 2>&1
mv $out.step1.fa $out.step1a.fa
ln -s $out.step1b.fa $out.step1.fa

";
$COMMAND="$COMMAND\n#END of WTDBG on $reads. Contigs can be found in $out.step1.fa\n\n
echo;echo PRIMARY CONTIG STATS:
perl $bin/seq_n50.pl $out.step1.fa
echo;date;echo FINISHED CSA STEP1;echo
";
print "$COMMAND\n";
