{
if(substr($1,1,7)=="ragout-")	{scf++}
else if($0!="" && substr($0,1,10)!="----------" && $1!="contig_1")	{
				i++;
				if($3==0 && $4!=""){$3=1};
				data[i]="scf"scf"\t"substr($1,2)"\t"substr($1,1,1)1"\t"$3;
				if($3==0 && $4=="" && guess!="on"){scf++};
				i++;
				data[i]="scf"scf"\t"substr($2,2)"\t"substr($2,1,1)1"\t"
				}
}

END{for(x=1;x<=i;x++)   {
			split(data[x],d,"\t");
			split(data[x+1],e,"\t");
			if(d[4]!="" || d[1]!=e[1])	{
							if(d[4]==""){d[4]=0;}
							#setting maxgap to 50000
							else if(d[4]>50000){d[4]=50000;}
							#setting min gap to 500
							else if(d[4]<500 && d[4]>=0 && d[1]==e[1]){d[4]=500;};
							print "-\t"d[1]"\t-\t"d[2]"\t-\t"d[3]"\t"d[4];
							}
			}
}
