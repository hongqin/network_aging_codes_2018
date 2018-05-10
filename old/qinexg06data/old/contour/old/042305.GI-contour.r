library(ade4);

tb <- read.table("new.new.merged.out",header=T);

expt.name = as.character( tb$expt );
strain = expt.name;
for( i in 1:length(expt.name)) {
   strain[i] = unlist( strsplit( expt.name[i], "\\.") ) [1];
}
strain;

xlim = c( 3E-4, max(tb$a) * 1.1 );
ylim = c( 0.05, max(tb$b) * 1.2 );
 
#plot( tb$b ~ tb$a ,  log='xy', col="red",  xlim=xlim, ylim=ylim );
plot( tb$b ~ tb$a ,  col="red",  xlim=xlim, ylim=ylim, xlab="I0",ylab="G" );
scatterutil.eti( tb$a, tb$b, strain, 0.5); 
 
 
  
 # I, G, mean-life-span, 2d contour plot, OK this is Gompertian Space
 
 predict.A <- function ( G, t ) {  G * log(0.5) / ( 1 - exp(G * t))  }
 
 step = 0.001; 
 G0 <- seq( 0.05, 0.3, by= step );   # i 
 t <- c( 25, 30, 35, 40, 45 );

 ii  <- matrix(nrow= length(G0), ncol= length(t) );
 
 for (j in 1:length(t) ) {
   for ( i in 1:length(G0) ) {
     ii[i,j] =  predict.A( G0[i], t[j] );
   }
 }

plot( tb$b ~ tb$a ,  col="red",  xlim=xlim, ylim=ylim, xlab="I0",ylab="G" );
for( i in 1:length(t)) { points ( ii[,i], G0, col="green", lty=2, type="l"); }
scatterutil.eti( tb$a, tb$b, strain, 0.5); 
title( main= "Natural replicative lifespan in Gompertzian space");

 