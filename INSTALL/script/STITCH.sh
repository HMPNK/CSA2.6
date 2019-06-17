export PATH=$2/bin:$PATH
$2/bin/seqtk seq -l 0 $1 | perl $2/script/last_scaffold_stitch.pl /dev/stdin| $2/bin/seqtk seq -l 100 -
