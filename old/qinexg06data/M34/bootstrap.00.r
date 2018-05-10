#092304Thu

## functions
 my.bootstrap <- function ( in.vector ) {
   data <- in.vector[ ! is.na(in.vector) ]; 
  ret = sample( data, size=length(data), replace=T);
 }
 
## the main section 
 tb <- read.table( "rls.tab");
 
 data <- tb[,1];
 data = data[ ! is.na(data) ]; 
 boot.out <- data.frame( matrix(nrow=length(data), ncol=length(data) ) );

 #go over each column and call my.bootstrap
 
 for ( j in 1:length(data) ) { 
  boot.out[,j] = my.bootstrap( data );
 }
 
 write.table( boot.out, "boot.tab", row.names=F, col.names=F, quote=F, sep="\t");
 

#q("no")

