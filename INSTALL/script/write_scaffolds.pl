#!/usr/bin/env perl
use strict;
use warnings;

sub usage
  {
    "usage: $0 contig-fasta-file scaffold-table\n";
  }

my ($fasta,$scftable)=@ARGV or die usage;

my %seq;
my $seqrc;
my $id;     
my @tab;
my@name;
my $oldscf="";

open(my $fastahandle, '<:encoding(UTF-8)', $fasta)
     or die "Could not open file '$fasta' $!";


#read from fasta file into hash:
while(my $row = <$fastahandle>) {
			chomp $row;
			if(substr($row,0,1) eq ">") {@name=split(/\s/,$row);$id=substr($name[0],1);$seq{$id}="";}
			else {$seq{$id}="$seq{$id}$row";}
			}

#read scaffold table and output scaffolds
open(my $scfhandle, '<:encoding(UTF-8)', $scftable)
     or die "Could not open file '$scftable' $!";

while(my $row = <$scfhandle>) {
				chomp $row;
				if($row ne "" && substr($row,0,1) ne "<" )	
					{
						@tab = split("\t",$row);
																			
						if($tab[1] ne $oldscf) {if($oldscf ne ""){print "\n"};print ">$tab[1]\n";};
						if($tab[5]==1) 
									{print "$seq{$tab[3]}"; } 
						
						else 			{
									$seqrc=reverse($seq{$tab[3]}) ;
									$seqrc=~tr/ACGTacgt/TGCAtgca/;
									print "$seqrc";
									$seqrc=""; 
									}
						print "N" x $tab[6];
						$oldscf=$tab[1];
					}
				}

print "\n";
