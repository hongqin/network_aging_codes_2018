# Gomertz plot for publication
#library(ade4);

nat <- read.table("nat.rls.tab", header=T);
lab <- read.table("lab.rls.tab", header=T);

isolates = as.character( nat$strain );
strains  = as.character( lab$strain );

#xlim = c( 2E-4, max(nat$a) * 1.1 );
#ylim = c( 0.05, max(nat$b) * 1.5 );
xlim = c( 2E-4, 0.02 );
ylim = c( 0.05, 0.25 );

x = nat$a;	y = nat$b - 0.005;
y[4]  = y[4] + 0.005 * 2;
y[14] = y[13] + 0.005 * 2;
x[14] = x[14] * 1.1;
x[13] = x[13] * 0.9;
 
postscript("gompertz-space.061405.ps", width=8, height=8)
plot( nat$b ~ nat$a ,  col="black",  xlim=xlim, ylim=ylim, 	xlab="I0",ylab="G", pch= 16, log='xy' );
text( x, y, isolates, cex=0.8); 

points( lab$b ~ lab$a ,  col="black", pch=1);
xx = lab$a; yy=lab$b - 0.005;
text( xx, yy, strains, cex=0.8 );
  
# ARLS contour lines
 predict.A <- function ( G, t ) {  G * log(0.5) / ( 1 - exp(G * t))  }
 
 step = 0.001; 
 G0 <- seq( 0.0525, 0.5, by= step );   # i 
 t <- c( 10, 15, 20, 25, 30, 35, 40, 45, 50, 60 );
 bb = 0.051
 aa <- ( log(0.5) * bb ) / (1 - exp(bb * t));
 ii  <- matrix(nrow= length(G0), ncol= length(t) );
 text ( aa, bb, t, cex=0.8); 

 for (j in 1:length(t) ) {
   for ( i in 1:length(G0) ) {
     ii[i,j] =  predict.A( G0[i], t[j] );
   }
 }

for( i in 1:length(t)) { points ( ii[,i], G0, col="black", lty=2, type="l"); }

#arrows( (tb$a - tb$a.sd), tb$b, (tb$a+tb$a.sd), tb$b, lwd=0.1, 	angle=90,code=3,lty=8, length=0.05,col="red");
#arrows( tb$a, (tb$b-tb$b.sd), tb$a, (tb$b+tb$b.sd), lwd=0.1, 	angle=90,code=3,lty=8, length=0.05,col="red");

dev.off();

 
