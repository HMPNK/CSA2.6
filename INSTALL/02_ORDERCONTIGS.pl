#!/usr/bin/env perl

#CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)

#THIS SCRIPT RUNS WHOLE GENOME ALIGNMENTS AND ORDERS CONTIGS

#INPUT DATA: 	Assembled contigs from 01_WTDBG or from other source 
#		suitable diverged REFERENCE GENOME OR backbone assembly


  use strict;
  use Getopt::Std;
my $userName =  $ENV{'LOGNAME'};
my $dir ="/home/$userName/CSA2.6";
my $script = "$dir/script";
my $bin = "$dir/bin";
my $contigs = "";
my $out = "CSA";
my $ref = "";
my $threads = 1;
my $guess= "off";

#______________________________________________________________________
#get options:
my %options=();
getopts("o:g:t:c:x:", \%options);

if(! $ARGV[0] && ! $options{c})       {
       print STDERR "
CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)
Author: Heiner Kuhl, Phd (kuhl\@igb-berlin.de)

THIS SCRIPT RUNS: assembly comparisons and contig ordering

INPUT DATA: Assembled contig, Reference sequences


      Options:\n
	      -c assembled sequences (as fasta file !)	
              -g reference sequences (suitable backbone assembly or reference genomes, as comma separated list)
              -o output-prefix
              -t number threads (default=1)
              -x on, off use guessed edges from RAGOUT (default=off), may improve chromosomal assembly on the cost of more structural errors
";       };
die "ERROR: NEED a contig file!!!\n" if(! $options{c});
die "ERROR: NEED at least one REFERENCE GENOME OR backbone for comparison!!!\n" if(! $options{g});
die "ERROR: unknown option $ARGV[0]" if $ARGV[0];

$contigs = $options{c} if($options{c});
$out = $options{o} if($options{o});
$ref = $options{g} if($options{g});
$threads = $options{t} if($options{t});
$guess = $options{x} if($options{x});

die "\nERROR: Invalid (-c) assembly file $contigs ...Exiting.\n\n" if(! -e $contigs);
die "\nERROR: Invalid (-c) suffix does not match .fa !!! Make sure to provide reads in gzipped fasta format! ...Exiting.\n\n" if( (substr($contigs, -3) ne ".fa") );
my @refcnt = split(/,/ , $ref);
foreach my $n (@refcnt) { die "$n does not exist ...Exiting." if(!-e $n);};

#02_CTGORDER:
my $COMMAND="#!/usr/bash\nset -e\nset -o pipefail\n\n";

$COMMAND="${COMMAND}\n#RUN CTGORDER\n\n";

my $maflist;
my $c;
foreach my $n (@refcnt) {
			my @m = split(/\//, $n);
			my $queryname = $m[-1];
			$queryname =~ s/\./_/g;
			$c++;

$COMMAND="${COMMAND}

echo;date;echo MASKING REFERENCE BY SELF-ALIGNMENT;echo

$bin/minimap2 -X -x asm5 -t $threads $n $n > $m[-1].SELF.paf 2>self-minimap.log
awk '{n=split(\$16,d,\":\");if(d[n]<=0.025){print \$1\"\\t\"\$3\"\\t\"\$4\"\\n\"\$6\"\\t\"\$8\"\\t\"\$9}}' $m[-1].SELF.paf | sort --buffer-size=16G --temporary-directory=./ -k1,1V -k2,2n | $bin/bedtools merge -d 100 -i - > $m[-1].mask.bed
$bin/seqtk seq $n | $bin/bedtools maskfasta -fi /dev/stdin -bed $m[-1].mask.bed -fo $m[-1].mask.fa

echo;date;echo CREATE LAST DATABASE FOR WHOLE GENOME ALIGNMENT;echo

$bin/lastdb -P $threads $m[-1]  $m[-1].mask.fa

echo;date;echo RUN WHOLE GENOME ALIGNMENT BY LAST;echo

bash $script/FASTLAST.sh $contigs $m[-1] $m[-1].maf $out $queryname $threads 1000000 $dir >> parallel.log 2>&1

echo;date;echo RUN RAGOUT TO ORDER CONTIGS ACCORDING TO REFERENCE;echo

awk 'BEGIN{print \".tree = ($out:0.01,$queryname:0.01);\\n.target = $out\\n.maf = $m[-1].maf\\n.blocks = 160000,80000,40000,20000,10000,5000\\n\\n$queryname.draft = true\\n$out.fasta = $contigs\"}' > $m[-1].recipe.txt

$bin/RAGOUT_V1.0/ragout.py  --no-refine --overwrite -o ./RAGOUT_$m[-1] -s maf $m[-1].recipe.txt >> ragout.log 2>&1

#CREATE NEW SCAFFOLDS DO SOME RENAMING ETC.
awk -v guess=$guess -f $script/ragout_link_2_tsv.awk ./RAGOUT_$m[-1]/scaffolds.links | awk '{gsub(\"scf\",\"R$c\_\",\$2);gsub(\" \",\"\\t\");print \$0}' > $out.scf$c.txt
$bin/seqtk seq -l 0 $contigs | perl $script/write_scaffolds.pl /dev/stdin $out.scf$c.txt | $bin/seqtk seq -l 100 - > $out.scf$c.fa
cut -f 4 $out.scf$c.txt > used$c.list
grep \">\" $contigs | grep -vw -F -f used$c.list  > unused$c.list
$bin/seqtk seq -l 0 $contigs | grep -A 1 -w -F -f unused$c.list| grep -vw ^\\-\\-\\\$  |$bin/seqtk seq -l 100 - >> $out.scf$c.fa
rm used$c.list unused$c.list			
";
			$contigs="$out.scf$c.fa";			
			$COMMAND="${COMMAND}\n";
			}
$COMMAND="${COMMAND}ln -s $out.scf$c.fa $out.step2.fa\n\n#END of CTGORDER\n
echo;echo PRIMARY SCAFFOLD STATS:
perl $bin/seq_n50.pl $out.step2.fa
echo;date;echo FINISHED CSA STEP2;echo
";
print "\n$COMMAND";
