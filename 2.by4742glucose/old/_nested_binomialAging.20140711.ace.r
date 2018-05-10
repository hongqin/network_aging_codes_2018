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
legend( 15.5,0.046, c("H0: t01 = t02,  loglikhood=-722.4", "H1: t01 != t02, loglikhood=-665.9") ) ;
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
legend( 17,0.05, c("H0 n1=n2 loglikhood=-700.52", "H1 n1!=n2 loglikhood=-665.9"), pch=c(20,20),col=c("black","gray")) ;
dev.off()



q("no")
#####

fitBinom1 = optim ( c(0.01,15, 3.5),  llh.binomialMortality.single.run , lifespan=x1, 
              method="L-BFGS-B", lower=c(1E-10,0.05,1E-1), upper=c(10,100,50) );  

fitBinom2 = optim ( c(0.01,10, 3.5),  llh.binomialMortality.single.run , lifespan=x2, 
                    method="L-BFGS-B", lower=c(1E-10,0.05,1E-1), upper=c(10,100,50) );  




H0( c(0.01, 15, 4, 0.01, 10, 4) )



#################################
##### LRT to exam whether two data on I and G.
### 1) model H0, same G and same I
### 2) model H1i, same G, different I
### 3) model H1g, different G, same I
### 4) model H2, different G, different I
                                  #    I1       G1       I2         G2
 H0  <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2] ) }  #all the same
 H1i <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[2] ) }  #different I
 H1g <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[4] ) }  # different G
 H2  <- function( rawIG ) { IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[4] ) }  # all different

 Hx.llh.gompertz.single.run <- function( rawIG, model, lifespan1, lifespan2 ) {
   IG = model(rawIG); 
   I1 = IG[1]; G1 = IG[2]; I2 = IG[3]; G2 = IG[4];
   my.data1 = lifespan1[!is.na(lifespan1)];
   my.data2 = lifespan2[!is.na(lifespan2)];
   log_s1 = (I1/G1) *(1 - exp(G1* my.data1))
   log_s2 = (I2/G2) *(1 - exp(G2* my.data2))
   log_m1 = log(I1) +  G1 * my.data1 ; 
   log_m2 = log(I2) +  G2 * my.data2 ; 
   my.lh = sum(log_s1) + sum(log_m1) + sum(log_s2) + sum(log_m2)
   print (IG ); #trace the convergence
   ret = - my.lh # because optim seems to maximize 
 }

## LRT to exam whether BY4742, BY4741 share the same G, I 
llh.H0  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=by4741, lifespan2=by4742 );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=by4741, lifespan2=by4742 );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=by4741, lifespan2=by4742 );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=by4741, lifespan2=by4742 );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );

## LRT to exam whether BY4743, BY4742 share the same G, I 
llh.H0  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=by4742, lifespan2=by4743 );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=by4742, lifespan2=by4743 );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=by4742, lifespan2=by4743 );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=by4742, lifespan2=by4743 );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);

LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );


#################################
##### LRT to exam THREE data on I and G.
### 1) model H0, same G and same I
### 2) model H1i, same G, different I
### 3) model H1g, different G, same I

 H0  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[1], rawIG[2]) }  #all the same

 H1g3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[1], rawIG[6]) } 

 H1i3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[5], rawIG[2]) } 

 H2ig3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[1], rawIG[2], rawIG[5], rawIG[6]) } 

 H3i2ig3  <- function( rawIG ) {
           #    I1       G1       I2         G2      I3           G3
  IG <- c(rawIG[1], rawIG[2], rawIG[3], rawIG[2], rawIG[5], rawIG[6]) } 

 H6  <- function( rawIG ) {  IG <- rawIG } 


Hx3.llh.gompertz.single.run <- function( rawIG, model, lifespan1, lifespan2, lifespan3 ) {
   IG = model(rawIG); 
   I1 = IG[1]; G1 = IG[2]; I2 = IG[3]; G2 = IG[4]; I3=IG[5]; G3=IG[6];
   my.data1 = lifespan1[!is.na(lifespan1)];
   my.data2 = lifespan2[!is.na(lifespan2)];
   my.data3 = lifespan3[!is.na(lifespan3)];
   log_s1 = (I1/G1) *(1 - exp(G1* my.data1))
   log_s2 = (I2/G2) *(1 - exp(G2* my.data2))
   log_s3 = (I3/G3) *(1 - exp(G3* my.data3))
   log_m1 = log(I1) +  G1 * my.data1 ; 
   log_m2 = log(I2) +  G2 * my.data2 ; 
   log_m3 = log(I3) +  G3 * my.data3 ; 
   my.lh = sum(log_s1) + sum(log_m1) + sum(log_s2) + sum(log_m2) + sum(log_s3) + sum(log_m3)
   print (IG ); #trace the convergence
   ret = - my.lh # because optim seems to maximize 
 }

llh.H0  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H0,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H1g3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H1g3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H1i3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H1i3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H2ig3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H2ig3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H3i2ig3  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H3i2ig3,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );

llh.H6  = optim( c(0.008,0.06,0.008,0.06, 0.01, 0.08), Hx3.llh.gompertz.single.run, model=H6,  lifespan1=by4741, lifespan2=by4742, lifespan3=by4743 );


cbind(llh.H0$par, llh.H1g3$par, llh.H1i3$par, llh.H2ig3$par, llh.H3i2ig3$par, llh.H6$par);
LH=cbind(llh.H0$value, llh.H1g3$value, llh.H1i3$value, llh.H2ig3$value, llh.H3i2ig3$value, llh.H6$value);
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );
1 - pchisq( 2*deltaLH, df =3 );


#################################
###### End of LRT    ############
#################################

