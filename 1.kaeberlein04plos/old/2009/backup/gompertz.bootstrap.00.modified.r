#092304Thu
 library(nlme)
 library(survival)

## functions
 gompfit <- function( data ) { #data is a vector
   tmp.fit <- survfit(Surv( data ) )
   s <- c( 1, tmp.fit$surv );
   t <- c( 0, tmp.fit$time );
   fm <- gnls( s ~ exp( (I/G)* ( 1 - exp(G * t) )  ) ,
                  start = list( I = 0.001, G = 0.2 ),
                  );
   I.fit <- fm$coefficients[1];
   G.fit <- fm$coefficients[2];
   c( I.fit, G.fit);
 }

 bootstrap.gompfit <- function ( in.vector ) {
   data <- in.vector[ ! is.na(in.vector) ]; 
   boot.num <- length(data);
   boot.out <- data.frame(matrix(ncol=3,nrow=boot.num));
   names(boot.out) <- c("I","G","avg.rls");
   for ( i in 1:boot.num ) {
     # data.tmp <- BY4741;
     data.tmp <- sample( data, size=length(data), replace=T);
     gf.out <- gompfit(data.tmp);
     boot.out[i, ] <- c( gf.out, mean(data.tmp) );
   }
   ret = boot.out;
 }
 
## the main section 
 rls <- read.table( "rls.tab", sep="\t", header=T);
 strains <- names(rls);
 
 #buffers for results
 return.labels <- c("I", "se.I", "G", "se.G", "rls", "se.rls" );
 gomp.out <- data.frame( matrix(nrow=length(names(rls)), ncol=length(return.labels) ) );
 names(gomp.out) <- return.labels;
 row.names(gomp.out) <- names(rls);

 #go over each column anb call bootstrap.gompfit
 group.1 = c(1:20 );
 group.2 = c(21:30 );
 group.3 = c(31:length(strains) );
  
 for ( j in group.1 ) { 
   boot.array = bootstrap.gompfit( rls[,j] );
   gomp.out[j,c(1,3,5)] = mean(boot.array);
   gomp.out[j,c(2,4,6)] = sd(boot.array);
 }

 #gomp.out;

 for ( j in group.2 ) { 
   boot.array = bootstrap.gompfit( rls[,j] );
   gomp.out[j,c(1,3,5)] = mean(boot.array);
   gomp.out[j,c(2,4,6)] = sd(boot.array);
 }

   
 for ( j in group.3 ) { 
   boot.array = bootstrap.gompfit( rls[,j] );
   gomp.out[j,c(1,3,5)] = mean(boot.array);
   gomp.out[j,c(2,4,6)] = sd(boot.array);
 }

  
 
 write.table(gomp.out, "092204.gomp.bootstrap.merged.lab-strains.csv", row.names=T, quote=T, sep="\t");
 

#q("no")


