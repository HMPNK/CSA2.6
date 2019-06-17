#!/usr/bin/env perl
#Heiner Kuhl 20181303
#Tool to find overlaps of neighbouring contig ends in scaffolds and close the gaps
#20190616 new version using last aligner

use strict;
use warnings;

sub usage
  {
    "usage: $0 <contig-fasta-file> \n";
  }

my ($fasta)=@ARGV or die usage;

my @file=split(/\//,$fasta);
if($file[-1] eq "stdin" || $file[-1] eq "-" || $file[-1] eq "STDIN") {$file[-1]="scaffold_stitch";}

my %seq;
my $id;
my @name;
my @namelist;
my $count = 0;
my $count2=0;

#change some parameters here if necessary!
my $max=30000;    # bp of contig ends to consider for overlaps (choose smaller values for short reads assemblies, 30000 is for long read sequencing)
my $maxmis1=100;  # bp of non aligned contig end allowed cut off 1
my $maxmis2=1000; # bp of non aligned contig end allowed cut off 2
my $minext=200;   # bp minimum length of contig A exrending contig B and vice versa (removes joining fully contained contig overlaps)
my $mins=100; # Score cut off; higher values are more stringent, but less gaps will be closed
my $minid=60; #minimum identity of last results

my $tseq;


open(OVL, ">$file[-1].closed.overlaps") or die "Could not open file '$file[-1].closed.overlaps' $!";
open(OVLB, ">$file[-1].skipped.overlaps") or die "Could not open file '>file[-1].closed.overlaps' $!";

open(my $fastahandle, '<:encoding(UTF-8)', $fasta)
     or die "Could not open file '$fasta' $!";


#read from fasta file into hash:
while(my $row = <$fastahandle>) {
                        chomp $row;
                        if(substr($row,0,1) eq ">") {@name=split(/\s/,$row);$id=substr($name[0],1);$seq{$id}="";$namelist[$count]=$id;$count=$count+1;}
                        else {$seq{$id}="$seq{$id}$row";}
                        }



for(@namelist) {
my $scfname;
my @scfcontigs;
my @scfgaps;

$scfname=$_;
@scfcontigs=split(/[Nn]+/,$seq{$_});

#if scaffold is a single contig, just output it to outfile! Else do the following stuff!
if(@scfcontigs==1) 	{
			print ">$_\n$seq{$_}\n";
			}
else{

#add IUPAC base support later!
@scfgaps=split(/[AaCcGgTt]+/,$seq{$_});

$count=0;
$count2=0;



for (my $i=0; $i+1 < @scfcontigs; $i++) {
		 open(OUT1, ">$scfname\_$count.fa");
		 my $len1 =length $scfcontigs[$count];
		 $tseq=substr $scfcontigs[$count] , -$max;
		 print OUT1 ">$scfname\_$count\n$tseq\n";
		 close(OUT1);
		 $count++;
		 $count2=$count-1;
		 open(OUT1, ">$scfname\_$count.fa");
		 my $len2 =length $scfcontigs[$count];
                 $tseq=substr $scfcontigs[$count] , 0 , $max;
		 print OUT1 ">$scfname\_$count\n$tseq\n";
                 close(OUT1);
		 $tseq="";
		
		
		`lastdb $scfname\_$count.fa $scfname\_$count.fa`;
		my $blast = `lastal -m 1 -a 1 -f BlastTab $scfname\_$count.fa  $scfname\_$count2.fa | grep -v '^#' | awk '{\$12=\$12+0;gsub(\" \",\"\\t\");print \$0\"\t$len1\t$len2\"}' | sort -k12,12rn`;
		`rm $scfname\_$count.fa* $scfname\_$count2.fa`;

#screen for end to end overlaps on plus strands with some non-overlapping regions allowed (larger indels)
#choose the best scoring end to end overlaps into list; remove small overlaps (low score/e-value)

		 my @ovllist=split(/\n/,$blast);
		 	for(@ovllist) 	{
					my @hit=split(/\t/,$_);
					if($hit[6]<$hit[7] && $hit[8]<$hit[9] ) {if($hit[12]>$max) {$hit[6]=$hit[6]+$hit[12]-$max;$hit[7]=$hit[7]+$hit[12]-$max;};

 										#good end to end overlap criteria here:
										if($hit[2]>=$minid && $hit[11]>=$mins && (($hit[8]<=$maxmis1 && $hit[12]-$hit[7]<=$maxmis2) || ($hit[8]<=$maxmis2 && $hit[12]-$hit[7]<=$maxmis1)) && $hit[6]>$hit[8]+$minext && $hit[13]-$hit[9]>$hit[12]-$hit[7]+$minext) 
											{
											for(@hit){print OVL "$_\t"};print OVL "\n";
											$scfcontigs[$count]=substr($scfcontigs[$count], $hit[9]);
										 	$scfcontigs[$count2]=substr($scfcontigs[$count2], 0, length($scfcontigs[$count2])-($hit[12]-$hit[7]));
											$scfgaps[$count]="";
											last;
											}
										else{for(@hit){print OVLB "$_\t"};print OVLB "\n";}
										 
										 
					}


		
				}        
#print scaffold again with trimmed contigs and removed gaps:
			}
	print ">$_\n";
	for (my $i=0; $i < @scfcontigs; $i++) {my $i2=$i; $i2++; if(! $scfgaps[$i2]){$scfgaps[$i2]=""};print "$scfcontigs[$i]$scfgaps[$i2]";}
	print "\n";

  
 }
}
