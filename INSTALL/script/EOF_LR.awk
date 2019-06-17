BEGIN{if(readl=="") {readl=500;}}
{
if(substr($1,1,1)==">"){
			n=$1
			} 

else if(length($1)>2*readl) {
				print n"\n"substr($1,1,readl);
				print n"\n"substr($1,length($1)-readl);
			}
}
