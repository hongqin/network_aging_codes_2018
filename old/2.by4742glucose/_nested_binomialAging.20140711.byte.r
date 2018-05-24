#test binomial aging model
rm(list=ls());

source("lifespan.20140711.r")

#################################
##### LRT to exam two datasets on binomial aging model, 
### model H0,  R, t0, n are all the same
### model H1i  R1<>R2
### model H1g, 
### model H2, 
#  rawIt0n = c( R1, t01, n1, R2, t02, n2)
H0  <- function( rawIt0n )  { Rt0n <- c(rawIt0n[1],rawIt0n[2],rawIt0n[3],rawIt0n[1],rawIt0n[2],rawIt0n[3]) }  #all the same
H1r  <- function( rawIt0n ) { Rt0n <- c(rawIt0n[1],rawIt0n[2],rawIt0n[3],rawIt0n[4],rawIt0n[2],rawIt0n[3]) }  # R1<>R2
H1t  <- function( rawIt0n ) { Rt0n <- c(rawIt0n[1],rawIt0n[2],rawIt0n[3],rawIt0n[1],rawIt0n[5],rawIt0n[3]) }  # t01 <> t02
H1n  <- function( rawIt0n ) { Rt0n <- c(rawIt0n[1],rawIt0n[2],rawIt0n[3],rawIt0n[1],rawIt0n[2],rawIt0n[6]) }  # n1 <> n2

Hx.llh.binomialMortality.single.run <- function( INrawRt0n, model, lifespan1, lifespan2 ) {
  Rt0n = model( INrawRt0n )
  I1 = Rt0n[1]; t01 = Rt0n[2]; n1=Rt0n[3];
  I2 = Rt0n[4]; t02 = Rt0n[5]; n2=Rt0n[6];  
  x1 = lifespan1[!is.na(lifespan1)];
  x2 = lifespan2[!is.na(lifespan2)];

  log_s1 = (I1 * t01 /n1 )*(1 - (1 + x1/t01)^n1);
  log_m1 = log(I1) +  (n1-1) * log(1 + x1/t01 ); 
  log_s2 = (I2 * t02 /n2 )*(1 - (1 + x2/t02)^n2);
  log_m2 = log(I2) +  (n2-1) * log(1 + x2/t02 ); 
  my.lh = sum(log_s1)  + sum(log_m1) + sum(log_s2) + sum(log_m2)
  print (Rt0n ); #trace the convergence
  ret = - my.lh # because optim seems to maximize 
}

set.seed(20140711)
x1 = round( rbinomial.aging( c(0.01, 15, 4), 100) + 0.5, digits=0);
x2 = round( rbinomial.aging( c(0.01, 10, 4), 100) + 0.5, digits=0);
llh.H0  = optim( c(0.01, 15, 4, 0.01, 15, 4), Hx.llh.binomialMortality.single.run, model=H0,lifespan1=x1, lifespan2=x2, 
                 method="L-BFGS-B", lower=c(1E-5,0.05,0.1, 1E-5,0.05,0.1), upper=c(10,100,50,10,100,50) );  
llh.H1t  = optim( c(0.01, 15, 4, 0.01, 10, 4), Hx.llh.binomialMortality.single.run, model=H1t,lifespan1=x1, lifespan2=x2, 
                 method="L-BFGS-B", lower=c(1E-5,0.05,0.1, 1E-5,0.05,0.1), upper=c(10,100,50,10,100,50) );  

cbind(llh.H0$par, llh.H1t$par)
LH = c(llh.H0$value, llh.H1t$value)
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );

#change of t0
set.seed(20140711)
x1 = round( rbinomial.aging( c(0.01, 25, 4.5), 100) + 0.5, digits=0); summary(x1)
x2 = round( rbinomial.aging( c(0.01, 24.0, 4.5), 100) + 0.5, digits=0); summary(x2)
llh.H0  = optim( c(0.01, 25, 4.5, 0.01, 25, 4.5), Hx.llh.binomialMortality.single.run, model=H0,lifespan1=x1, lifespan2=x2, 
                 method="L-BFGS-B", lower=c(1E-5,0.05,0.1, 1E-5,0.05,0.1), upper=c(10,100,50,10,100,50) );  
llh.H1n  = optim( c(0.01, 25, 4.5, 0.01, 22, 4.5), Hx.llh.binomialMortality.single.run, model=H1n,lifespan1=x1, lifespan2=x2, 
                  method="L-BFGS-B", lower=c(1E-5,0.05,0.1, 1E-5,0.05,0.1), upper=c(10,100,50,10,100,50) );  
cbind(llh.H0$par, llh.H1n$par)
LH = c(llh.H0$value, llh.H1t$value); LH
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );

my.breaks = seq( 0,50,by=5);
h1 <- hist(x1, br= my.breaks, xlab = "replicative lifespan", ylab = "relative density", freq=F ) ;
h2 <- hist(x2, br= my.breaks, xlab = "replicative lifespan", ylab = "relative density", freq=F ) ;
x12 <-  data.frame( rbind(h1$density,h2$density) )  ;
names( x12 ) <- h1$mids;
row.names(x12) <- c( "R=0.01 t0=25 n=4.5", "R=0.01 t0=24 n=4.5" )
pdf("_LRT-R0.01-t0-25vs24-n4.5--20140711.pdf")
barplot( as.matrix(x12), beside=T, col=c("black","gray"), ylab="Relative Density", xlab="Replicative Lifespan")
legend( 17,0.045, c("t0=25.0 loglikhood=-722.4", "t0=24.0 loglikhood=-665.9"), pch=c(20,20),col=c("black","gray")) ;
text(21.5, 0.037, "LRT p<0.001")
dev.off()



#change of n
set.seed(20140711)
x1 = round( rbinomial.aging( c(0.01, 25, 6), 100) + 0.5, digits=0); summary(x1)
x2 = round( rbinomial.aging( c(0.01, 25, 4.5), 100) + 0.5, digits=0); summary(x2)
llh.H0  = optim( c(0.01, 25, 6, 0.01, 25, 6), Hx.llh.binomialMortality.single.run, model=H0,lifespan1=x1, lifespan2=x2, 
                 method="L-BFGS-B", lower=c(1E-5,0.05,0.1, 1E-5,0.05,0.1), upper=c(10,100,50,10,100,50) );  
llh.H1n  = optim( c(0.01, 25, 6, 0.01, 25, 3), Hx.llh.binomialMortality.single.run, model=H1n,lifespan1=x1, lifespan2=x2, 
                  method="L-BFGS-B", lower=c(1E-5,0.05,0.1, 1E-5,0.05,0.1), upper=c(10,100,50,10,100,50) );  

cbind(llh.H0$par, llh.H1n$par)
LH = c(llh.H0$value, llh.H1t$value); LH
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );

my.breaks = seq( 0,50,by=5);
h1 <- hist(x1, br= my.breaks, xlab = "replicative lifespan", ylab = "relative density", freq=F ) ;
h2 <- hist(x2, br= my.breaks, xlab = "replicative lifespan", ylab = "relative density", freq=F ) ;
x12 <-  data.frame( rbind(h1$density,h2$density) )  ;
names( x12 ) <- h1$mids;
row.names(x12) <- c( "R=0.01 t0=25 n=6", "R=0.01 t0=25 n=4.5" )
pdf("_LRT-R0.01-t0-25-n6-4.5-20140711.pdf")
barplot( as.matrix(x12), beside=T, col=c("black","gray"), ylab="Relative Density", xlab="Replicative Lifespan")

legend( 19,0.05, c("n=6 loglikhood=-700.52", "n=4.5 loglikhood=-665.9"), pch=c(20,20),col=c("black","gray")) ;
text(22, 0.042, "LRT p<0.001")
legend( 17,0.05, c("H0 n1=n2 loglikhood=-700.52", "H1 n1!=n2 loglikhood=-665.9"), pch=c(20,20),col=c("black","gray")) ;
dev.off()

