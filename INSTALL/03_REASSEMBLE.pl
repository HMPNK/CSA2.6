#!/usr/bin/env perl

#CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)

#THIS SCRIPT RUNS extracts reads around gaps and reassembles them locally, finally it reassembles the whole genome with the potentially gap closing contigs

#INPUT DATA: pacbio or oxford nanopore long reads (provided as fa.gz) and scaffolded genome assembly


  use strict;
  use Getopt::Std;
my $userName =  $ENV{'LOGNAME'};
my $dir ="/home/$userName/CSA2.6";
my $wtdbg = "$dir/bin/wtdbg2.2";
my $bin = "$dir/bin";
my $script = "$dir/script";
my $reads = "";
my $out = "CSA";
my $threads = 1;
my $assemblypath ;

#______________________________________________________________________
#get options:
my %options=();
getopts("o:r:a:t:", \%options);

if(! $ARGV[0] && ! $options{r})       {
       print STDERR "
CSA version 2.6c (Chromosome Scale Assembler, employing synteny with diverged reference genomes)
Author: Heiner Kuhl, Phd (kuhl\@igb-berlin.de)

THIS SCRIPT extracts reads mapped near gaps and runs WTDBG gap re-assembly and assembly of gap contigs with the input assembly


INPUT DATA: pacbio or oxford nanopore long reads (provided as fa.gz)


      Options:\n
              -a scaffolded wtdbg assembly
              -r input long reads as fa.gz (required!)
              -o output-prefix
              -t number threads (default=1) 

";       };
die "ERROR: NEED a read file!!!\n" if(! $options{r});
die "ERROR: NEED a scaffolded assembly file!!!\n" if(! $options{a});
die "ERROR: unknown option $ARGV[0]" if $ARGV[0];

$reads = $options{r} if($options{r});
$out = $options{o} if($options{o});
$assemblypath = $options{a} if($options{a});
my @m = split(/\//, $assemblypath);
                        my $assembly = $m[-1];


$threads = $options{t} if($options{t});


die "\nERROR: Invalid (-r) long read file $reads ...Exiting.\n\n" if(! -e $reads);
die "\nERROR: Invalid (-r) does not match '.fa.gz' !!! Make sure to provide reads in gzipped fasta format! ...Exiting.\n\n" if(substr($reads, -5) ne "fa.gz");



#03_REASSEMBLEGAPS:
my $COMMAND="#!/usr/bash\nset -e\nset -o pipefail\n\n";

$COMMAND="$COMMAND\n#RUN GAP REASSEMBLY by $reads and $assembly\n\n";
$COMMAND="$COMMAND
#create bed file for Gaps in the assembly
ln -s $assemblypath $assembly
$bin/samtools faidx $assembly
cut -f1,2 $assembly.fai  > $assembly.sizes
$bin/hgFakeAgp -minContigGap=0 -minScaffoldGap=0 $assembly $assembly.agp
grep -w N $assembly.agp |awk '{print \$1\"\\t\"\$2-1\"\\t\"\$3}'  > $assembly.gap.bed
grep -w D $assembly.agp |awk '{print \$1\"\\t\"\$2-1\"\\t\"\$3}'  > $assembly.contigs.bed

# add 20000bp to each side of gap
awk -v flank=20000 '{if(\$2-flank<0) {\$2=0} else {\$2=\$2-flank};\$3=\$3+flank;print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$1\"-gap\"x[\$1]++}' $assembly.gap.bed > $assembly.gap_flanked20kbp.bed
# add 20000bp of each scaffold end and non scafollded contig end
awk -v flank=20000 '{if(\$2>2*flank){print \$1\"\\t0\\t\"flank\"\\t\"\$1\"-start\";print \$1\"\\t\"\$2-flank\"\\t\"\$2\"\\t\"\$1\"-end\";} else {print \$1\"\\t0\\t\"\$2\"\\t\"\$1\"-start\";}}' $assembly.sizes >> $assembly.gap_flanked20kbp.bed

echo;date;echo MAP RAW READS TO ASSEMBLY;echo

#USE MINIMAP2 PAF mode, it is faster and uses much less CPU than KBM2!
$bin/minimap2 -I 100G -t $threads -x map-pb $assembly $reads 2>$assembly.minimap2.log > $assembly.paf

#create bed with read positions
awk '{if(\$12>=20 && o!=\$1){print \$6\"\\t\"\$8\"\\t\"\$9\"\\t\"\$1};o=\$1}' $assembly.paf > $assembly.read_placement.bed

#intersect with flanked gaps
$bin/bedtools intersect -wo -b $assembly.gap_flanked20kbp.bed -a $assembly.read_placement.bed | awk '{print \$4\"\\t\"\$8}' > $assembly.reads_to_gaps.list

echo;date;echo WRITE GAP/CTGEND MATCHING RAW READS TO FILES;echo

$bin/pigz -dc $reads | $bin/seqtk seq -l 0 - | awk -v infile=$assembly.reads_to_gaps.list 'BEGIN{while(getline l < infile){split(l,d,\"\\t\");d[1]=\">\"d[1];gap[d[1]]=d[2]}} {if(substr(\$1,1,1)==\">\"){outfile=gap[\$1]\".fasta\"};if(outfile!=\".fasta\"){print \$0 >> outfile};close(outfile)}'

echo;date;echo RUN LOCAL ASSEMBLIES;echo


cut -f 2 $assembly.reads_to_gaps.list > greplist
#use xargs here to not run into file number limits!!!
ls | grep -F -f greplist | grep fasta | xargs wc -l | sort -k1,1rn |awk '{i++;if(\$1>=4 && i>1){print \"sh $script/wtdbg-BATCH.sh \"\$2\" $wtdbg\/\"}}' > $assembly.wtdbg-BATCH-ALL-RUNS.sh
#run local reassemblies in parallel, timeout is now implemented in wtdbg-BATCH.sh -> sometimes happens with problematic repeats
cat $assembly.wtdbg-BATCH-ALL-RUNS.sh| $bin/parallel -j $threads > parallel.log 2>&1

#merge contigs
ls | grep -F -f greplist | grep wtdbg.fa | xargs cat > $assembly.GAP_HELPERS.fa
#make 5X coverage for GAP_Helpers (and add  unique fasta ids!)
cat $assembly.GAP_HELPERS.fa $assembly.GAP_HELPERS.fa $assembly.GAP_HELPERS.fa $assembly.GAP_HELPERS.fa $assembly.GAP_HELPERS.fa| $bin/seqtk seq -l 0 - |awk '{if(substr(\$1,1,1)==\">\") {n=\$1\"_\"c[\$1]++} else {print n\"\\n\"\$1}}'| $bin/pigz -c > $assembly.GAP_HELPERS_5X.fa.gz
rm $assembly.GAP_HELPERS.fa
#write gap assembly logs to single file
ls | grep -F -f greplist | grep wtdbg.log | xargs cat> $assembly.GAP.assemblies.log
#remove single files
ls| grep -F -f greplist | xargs rm


#create Fake reads for assembly with gap fillers

$bin/bedtools makewindows -w 155555 -b $assembly.contigs.bed  > $assembly.fakeReads.bed
$bin/bedtools makewindows -w 177777 -b $assembly.contigs.bed >> $assembly.fakeReads.bed
$bin/bedtools makewindows -w 199999 -b $assembly.contigs.bed >> $assembly.fakeReads.bed
$bin/bedtools makewindows -w 222222 -b $assembly.contigs.bed >> $assembly.fakeReads.bed
$bin/bedtools makewindows -w 233333 -b $assembly.contigs.bed >> $assembly.fakeReads.bed
$bin/bedtools makewindows -w 255555 -b $assembly.contigs.bed >> $assembly.fakeReads.bed

#reverse complement
$bin/seqtk seq -r $assembly > $assembly.RC
$bin/hgFakeAgp -minContigGap=0 -minScaffoldGap=0 $assembly.RC $assembly.RC.agp
grep -w D $assembly.RC.agp |awk '{print \$1\"\\t\"\$2-1\"\\t\"\$3}'  > $assembly.RC.contigs.bed

$bin/bedtools makewindows -w 155555 -b $assembly.RC.contigs.bed  > $assembly.RC.fakeReads.bed
$bin/bedtools makewindows -w 177777 -b $assembly.RC.contigs.bed >> $assembly.RC.fakeReads.bed
$bin/bedtools makewindows -w 199999 -b $assembly.RC.contigs.bed >> $assembly.RC.fakeReads.bed
$bin/bedtools makewindows -w 222222 -b $assembly.RC.contigs.bed >> $assembly.RC.fakeReads.bed
$bin/bedtools makewindows -w 233333 -b $assembly.RC.contigs.bed >> $assembly.RC.fakeReads.bed
$bin/bedtools makewindows -w 255555 -b $assembly.RC.contigs.bed >> $assembly.RC.fakeReads.bed


$bin/bedtools getfasta -bed $assembly.fakeReads.bed -fi $assembly -fo $assembly.fakeReads_6X.fasta > /dev/null 2>&1
$bin/bedtools getfasta -bed $assembly.RC.fakeReads.bed -fi $assembly.RC -fo $assembly.fakeReads_6X.fasta.RC > /dev/null 2>&1
cat $assembly.fakeReads_6X.fasta $assembly.fakeReads_6X.fasta.RC |awk '{if(substr(\$1,1,1)==\">\"){print \">\"x++} else{print}}' | $bin/pigz -c  > $assembly.fakeReads_12X.fasta.gz
rm $assembly.fakeReads_6X.fasta $assembly.fakeReads_6X.fasta.RC

echo;date;echo REASSEMBLE WHOLE GENOME USING REASSEMBLED GAP/CTGENDS AND PRIOR CONTIGS;echo

$bin/pigz -dc $assembly.fakeReads_12X.fasta.gz $assembly.GAP_HELPERS_5X.fa.gz | $wtdbg/wtdbg2 -t $threads -i /dev/stdin -o $assembly.GAPFILL-WTDBG -p 23 -S 4 -s 0.5 -e 6 -A 2>$assembly.GAPFILL-wtdbg.log

#use fast mode consenser here, because fake reads have already computed by higher qualty consener! Use -S 0 otherwise crashed on very large contigs (Try to find better solution)
$wtdbg/wtdbg-cns -i $assembly.GAPFILL-WTDBG.ctg.lay.gz -S 0 -f -t $threads -o $assembly.GAPFILL.fa 2>>$assembly.GAPFILL-wtdbg.log

echo;date;echo DO SOME CHECKS AND REMOVE POTENTIAL MISSASSEMBLIES;echo

$bin/minimap2 -I 100G -t $threads -x asm5 $assembly.GAPFILL.fa $assembly > $assembly.GAPFILL.paf 2>>$assembly.minimap2.log
#INTERSCAFFOLD
awk '{if(\$12>=60 && \$9-\$8 > 30000){print \$1\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9}}' $assembly.GAPFILL.paf| sort -k2,2V -k3,3n -k4,4n| grep -v ^ctg|awk '{if(\$2==oq && \$1!=oR){n=split(ol,d,\"\\t\");if(d[n]<=\$3){print \$2\"\\t\"d[n]\"\\t\"\$3} else {print \$2\"\\t\"\$3\"\\t\"d[n]}};oR=\$1;oq=\$2;ol=\$0}' | awk '{if(\$2<=\$3){print;} else{print \$1\"\\t\"\$3\"\\t\"\$2;}}' > $assembly.GAPFILL.SPLIT.bed
#INTRASCAFFOLD
awk '{if(\$12>=60 && \$9-\$8 > 30000){if(\$5==\"+\"){print \$1\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9\"\\t\"\$3\"\\t\"\$4\"\\t\"\$5} else if(\$5==\"-\"){print \$1\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9\"\\t\"\$4\"\\t\"\$3\"\\t\"\$5}}}' $assembly.GAPFILL.paf|grep -v ^ctg|sort -k2,2V -k3,3n -k4,4n|awk '{if(o1==\$1 && o2==\$2){split(l,d,\"\\t\"); print \$1\"\\t\"\$2\"\\t\"d[4]\"\\t\"\$3\"\\t\"d[6]\"\\t\"\$5\"\\t\"\$7\"\\tGAP\\t\"\$3-d[4]\"\\t\"\$5-d[6]};print;o1=\$1;o2=\$2;l=\$0}'|awk '{if(\$10^2>60000^2 || \$9^2>60000^2 && \$0!=\"\"){if((\$9/\$10)^2>10^2 || (\$10/\$9)^2>10^2)print \$0}}'| cut -f 2-4 | awk '{if(\$2<=\$3){print;} else{print \$1\"\\t\"\$3\"\\t\"\$2;}}' >> $assembly.GAPFILL.SPLIT.bed
#CONTIG ENDS
awk '{if(\$12>=60 && \$9-\$8 > 30000){print \$1\"\\t\"\$6\"\\t\"\$8\"\\t\"\$9}}' $assembly.GAPFILL.paf | sort -k2,2V -k3,3n -k4,4n| awk '{i++;h[i]=\$0} END{for(x=1;x<=i;x++){n1=split(h[x-1],a,\"\\t\");n2=split(h[x],b,\"\\t\");n3=split(h[x+1],c,\"\\t\");if(substr(b[1],1,3)==\"ctg\" && b[2]==c[2] && a[2]!=b[2]) {print b[2]\"\\t\"b[4]\"\\t\"c[3]} else if(substr(b[1],1,3)==\"ctg\" && a[2]==b[2] && b[2]!=c[2]) {print b[2]\"\\t\"a[4]\"\\t\"b[3]}}}' | awk '{if(\$2<=\$3){print;} else{print \$1\"\\t\"\$3\"\\t\"\$2;}}' >> $assembly.GAPFILL.SPLIT.bed

sort -k1,1V -k2,2n -k3,3rn $assembly.GAPFILL.SPLIT.bed > temp
mv temp $assembly.GAPFILL.SPLIT.bed

if [ -s $assembly.GAPFILL.SPLIT.bed ]
then
$bin/samtools faidx $assembly.GAPFILL.fa
cut -f 1,2 $assembly.GAPFILL.fa.fai > $assembly.GAPFILL.sizes
$bin/bedtools complement -g $assembly.GAPFILL.sizes -i $assembly.GAPFILL.SPLIT.bed > $assembly.GAPFILL.GOOD.bed
else
$bin/samtools faidx $assembly.GAPFILL.fa
cut -f 1,2 $assembly.GAPFILL.fa.fai > $assembly.GAPFILL.sizes
awk '{print \$1\"\\t0\\t\"\$2}' $assembly.GAPFILL.sizes > $assembly.GAPFILL.GOOD.bed
fi

sort -k1,1V -k2,2n $assembly.GAPFILL.SPLIT.bed $assembly.GAPFILL.GOOD.bed | $bin/bedtools merge -d -1000 | awk '{if(\$3-\$2>=5000){print \$0\"\\t\"\$1\"_\"x[\$1]++}}' > $assembly.GAPFILL.REGIONS.bed
$bin/bedtools getfasta -name -fi $assembly.GAPFILL.fa -fo $out.step3.fa -bed $assembly.GAPFILL.REGIONS.bed > /dev/null 2>&1
";
$COMMAND="$COMMAND\n#END of GAP REASSEMBLY and CLOSURE. IMPROVED CONTIGS can be found in: $out.step3.fa\n\n
echo;echo STEP3 GAP CLOSED CONTIG STATS:
perl $bin/seq_n50.pl $out.step3.fa
echo;date;echo FINISHED CSA STEP3;echo

";
print "$COMMAND\n";
