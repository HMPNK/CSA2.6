##USAGE: bash FASTLAST.sh <query-fasta> <lastdb> <maf-outname> <queryspecies-id> <targetspecies-id> <cpu-threads> <splitted query size> <bin/script path>

##FOR very large query contigs LAST aligner will run slow as parallelization is not working efficiently!
##workaround
##create splitted query
ln -s $1 ./FASTLAST

$8/bin/seqtk comp FASTLAST | cut -f 1,2 > FASTLAST.sizes
$8/bin/bedtools makewindows -w $7 -g FASTLAST.sizes > FASTLAST.split.bed
$8/bin/bedtools getfasta -bed FASTLAST.split.bed -fi FASTLAST -fo FASTLAST.split.fa
##align
LAST="$8/bin/lastal -m 1 -l 1 -P 1 $2 - | $8/bin/last-split -m 0.0001 | $8/bin/maf-swap"
cat FASTLAST.split.fa \
| $8/bin/parallel -j $6 --recstart '>' --block-size $7 --pipe \
$LAST \
| $8/bin/last-split -m 0.0001 \
| grep -vw ^p \
|awk -f $8/script/prepare_maf.awk -v query=$5 -v target=$4 \
|awk -v chromsizes=FASTLAST.sizes -f $8/script/update_maf_coords.awk > $3

##remove intermediate fasta
rm FASTLAST*
