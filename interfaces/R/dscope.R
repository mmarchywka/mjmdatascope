
#Rcpp::compileAttributes("mjmdscope")
library("mjmdscope")
nr=100.0
nc=5
df<- as.data.frame(matrix(ncol=nc,nrow=nr,data=0.0))
#df<-data.frame(nrows=nr,cols=nc)
for(k in 1:10000)
{
for(i in 1.0:nr)
{

for(j in 1.0:(nc-1))
{
df[i,j]=sin(.01*2.0*pi*j*i/nc+.0333*k*2*pi);
}  # j
df[i,nc]=1.0*(i)/nr;
} # i 
mjmdscope::mjmdscopesend(df)
Sys.sleep(.2)
} # k 
str(df)
