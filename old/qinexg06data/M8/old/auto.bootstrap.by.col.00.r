#040405Mon

 my.bootstrap <- function ( in.vector ) {
   data <- in.vector[ ! is.na(in.vector) ]; 
   boot.num <- ( length(data) %/% 10 ) * 10;
   boot.out <- data.frame(matrix(ncol=boot.num, nrow=length(data) ));
   for ( j in 1:boot.num ) {
     data.tmp <- sample( data, size=length(data), replace=T);
     data.tmp2 <- sort( data.tmp );
     boot.out[ , j ] <- data.tmp2;
   }
   ret = boot.out;
 }
 
## the main section 
# rls <- read.table( "092104.naturals.rls.all.tab", sep="\t", header=T);
 rls <- read.table( "rls.tab", sep="\t", header=T);
 in.col.num = length( names(rls) );
 in.vector = unlist( rls );
 
 bootstrap = my.bootstrap( in.vector ); 

 write.table( bootstrap, "bootstrap.rls.tab", sep="\t", row.names=F, col.names=F);
 
#q("no")
