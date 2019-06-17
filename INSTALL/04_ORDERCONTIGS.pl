#!/usr/bin/env perl

#CSA version 2.6 (Chromosome Scale Assembler, employing synteny with diverged reference genomes)

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
my @m;

#______________________________________________________________________
#get options:
my %options=();
getopts("o:g:t:c:x:", \%options);

if(! $ARGV[0] && ! $options{c})       {
       print STDERR "
CSA version 2.6 (Chromosome Scale Assembler, employing synteny with diverged reference genomes)
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
my $contigs2 = $contigs;
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
			@m = split(/\//, $n);
			my $queryname = $m[-1];
			$queryname =~ s/\./_/g;
			$c++;

$COMMAND="${COMMAND}
#CREATE LAST DATABASE FOR WHOLE GENOME ALIGNMENT
ln -s ../02_ORDERCTG/$m[-1].??? .
rm $m[-1].maf

echo;date;echo RUN WHOLE GENOME ALIGNMENT BY LAST;echo

bash $script/FASTLAST.sh $contigs $m[-1] $m[-1].maf $out $queryname $threads 1000000 $dir

echo;date;echo RUN RAGOUT TO ORDER CONTIGS ACCORDING TO REFERENCE;echo

awk 'BEGIN{print \".tree = ($out:0.01,$queryname:0.01);\\n.target = $out\\n.maf = $m[-1].maf\\n.blocks = 80000,60000,40000,20000\\n\\n$queryname.draft = true\\n$out.fasta = $contigs\"}' > $m[-1].recipe.txt

$bin/RAGOUT_V1.0/ragout.py  --no-refine --overwrite -o ./RAGOUT_$m[-1] -s maf $m[-1].recipe.txt

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

$COMMAND="${COMMAND}
#sort fasta by size 
$bin/seqkit sort -lr $out.scf$c.fa > temp.fa
mv temp.fa $out.scf$c.fa

#RUN FINAL GAP CLOSE (find contig end overlaps by blast and join)
bash $script/STITCH.sh $out.scf$c.fa $dir > ../$out.final.A.fa
mv scaffold_stitch.closed.overlaps scaffold_stitch.closed.overlaps.A; mv scaffold_stitch.skipped.overlaps scaffold_stitch.skipped.overlaps.A
";
if($guess eq "off") {
if($c>1) { $c=$c-1; $contigs="$out.scf$c.fa"; } else { $contigs=$contigs2; };
$COMMAND="${COMMAND}
#CREATE ASSEMBLY B using RAGOUTS guessed edges
awk -v guess=on -f $script/ragout_link_2_tsv.awk ./RAGOUT_$m[-1]/scaffolds.links | awk '{gsub(\"scf\",\"R$c\_\",\$2);gsub(\" \",\"\\t\");print \$0}' > $out.scf$c.B.txt
$bin/seqtk seq -l 0 $contigs | perl $script/write_scaffolds.pl /dev/stdin $out.scf$c.B.txt | $bin/seqtk seq -l 100 - > $out.scf$c.B.fa
cut -f 4 $out.scf$c.B.txt > used$c.B.list
grep \">\" $contigs | grep -vw -F -f used$c.B.list  > unused$c.B.list
$bin/seqtk seq -l 0 $contigs | grep -A 1 -w -F -f unused$c.B.list| grep -vw ^\\-\\-\\\$  |$bin/seqtk seq -l 100 - >> $out.scf$c.B.fa
rm used$c.B.list unused$c.B.list

#sort fasta by size
$bin/seqkit sort -lr $out.scf$c.B.fa > temp.fa
mv temp.fa $out.scf$c.B.fa

echo;date;echo RUN FINAL GAP CLOSE (find contig end overlaps by blast and join);echo
bash $script/STITCH.sh $out.scf$c.B.fa $dir > ../$out.final.B.fa
mv scaffold_stitch.closed.overlaps scaffold_stitch.closed.overlaps.B; mv scaffold_stitch.skipped.overlaps scaffold_stitch.skipped.overlaps.B
";
}

$COMMAND="${COMMAND}
echo;date;echo FINISHED CSA STEP4, CSA COMPLETED!;echo
\n#END of CTGORDER\n
";


print "\n$COMMAND";
