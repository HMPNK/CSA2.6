BEGIN{

while(getline l < chromsizes)
{split(l,d,"\t");csize[d[1]]=d[2];}

print "#\tLAST\tversion 880";
}

{
if($1=="a"){print;i=0}
else if($1=="s" && i==0){
			n=split($2,a,":");
			split(a[n],b,"-");
			split(a[1],c,"\.");

			n=split($0,d,"\t");
			d[2]=a[1];
			d[6]=csize[c[2]];
			d[3]=d[3]+b[1];
			for(x=1;x<n;x++) {printf d[x]"\t"};printf d[n]"\n"
			i++;
			}

else if($1=="s" && i>0)	{
			print $0"\n";
			}
}
