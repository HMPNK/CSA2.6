#H.Kuhl, kuhl@igb-berlin.de

BEGIN{
if(target==0){target="target"}
if(query==0){query="query"}
}

{
n=split($0,d," ");

if(d[1]=="s") {
		i++;
		if(i==1) {d[2]=target"."d[2]};
		if(i==2) {d[2]=query"."d[2];i=0;};
		}

for(x=1;x<n;x++) {printf d[x]"\t";}
printf d[n]"\n";

}

