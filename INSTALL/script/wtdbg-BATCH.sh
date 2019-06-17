timeout 20m $2wtdbg2 -f -t 1 -i $1 -o $1-WTDBG -K 3000 -p 19 -S 2 -e 2 -A 2>$1-wtdbg.log
timeout 10m $2wtdbg-cns -f -i $1-WTDBG.ctg.lay.gz -t 1 -o $1-wtdbg.fa 2>>$1-wtdbg.log
rm $1-WTDBG*
