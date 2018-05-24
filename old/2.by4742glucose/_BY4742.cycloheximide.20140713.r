#test binomial aging model on cycloheximide dose in BY4742

rm(list=ls());
setwd("~/projects/0.network.aging.prj/0a.network.fitting")
source("lifespan.20140711.r")
require(flexsurv)

#files = list.files(, path='DR', pattern="cyclo")
mixed_names = c(
  "BY4742_MATalpha_temp30_media_YPD_n19930.csv", 0,
  "BY4742_MATalpha_temp30_media_YPD  10 ngmL cycloheximide_n80.csv",  10,
#  "BY4742_MATalpha_temp30_media_YPD + 10ngml cycloheximide_n120.csv", 10,
  "BY4742_MATalpha_temp30_media_YPD  10ngml cycloheximide_n120.csv",  10,
  "BY4742_MATalpha_temp30_media_YPD + 25 ngmL cycloheximide_n160.csv", 25,
  "BY4742_MATalpha_temp30_media_YPD  25 ngmL cycloheximide_n160.csv", 25,
#  "BY4742_MATalpha_temp30_media_YPD  30ngml cycloheximide_n40.csv",   30,
#  "BY4742_MATalpha_temp30_media_YPD  35ngml cycloheximide_n40.csv",   35,
#  "BY4742_MATalpha_temp30_media_YPD  40ngml cycloheximide_n40.csv",   40,
#  "BY4742_MATalpha_temp30_media_YPD  45ngml cycloheximide_n40.csv",   45,
  "BY4742_MATalpha_temp30_media_YPD  50 ngmL cycloheximide_n37.csv",  50,
  "BY4742_MATalpha_temp30_media_YPD  50ngml cycloheximide_n140.csv",  50,
#  "BY4742_MATalpha_temp30_media_YPD + 30ngml cycloheximide_n40.csv",  30,
#  "BY4742_MATalpha_temp30_media_YPD + 35ngml cycloheximide_n40.csv",  35,
#  "BY4742_MATalpha_temp30_media_YPD + 40ngml cycloheximide_n40.csv",  40,
#  "BY4742_MATalpha_temp30_media_YPD + 45ngml cycloheximide_n40.csv",  45,
  "BY4742_MATalpha_temp30_media_YPD + 50ngml cycloheximide_n140.csv", 50,
#  "BY4742_MATalpha_temp30_media_YPD + 100ngml cycloheximide_n139.csv",100,
  "BY4742_MATalpha_temp30_media_YPD  100ngml cycloheximide_n139.csv",100
);
files = mixed_names[seq(1,(length(mixed_names)-1), 2)]
cycloheximide = as.numeric( mixed_names[seq(2,(length(mixed_names)), 2)] )

output = data.frame(files); 
output$mean = NA; output$stddev = NA; output$CV=NA; 
output$R=NA; output$t0=NA; output$n=NA;
output$Rflex = NA; output$Gflex = NA; 
output$cycloheximide = cycloheximide

pdf("BY4742CyclohexmideDose-binomial-20140713.pdf", width=10, height=8)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))

#i=1; 
for( i in 1:length(files)) {
  filename = paste('DR/', files[i], sep='')
  tb = read.csv(filename)
  summary(tb)
  output$mean[i]=mean(tb$rls)
  output$stddev[i]= sd(tb$rls)
  output$CV[i] = sd(tb$rls) / mean(tb$rls)
  
  fitGomp = flexsurvreg( formula=Surv(tb$rls) ~ 1, dist='gompertz')
  Rhat = fitGomp$res[2,1]; Ghat = fitGomp$res[1,1] 
  output$Rflex[i] = Rhat; output$Gflex[i]=Ghat; 
  
  #G = (n-1) / t0, so t0 = (n-1)/G
  nhat = 3; t0= (nhat-1)/Ghat; t0
  fitBinom = optim ( c(Rhat, t0, nhat),  llh.binomialMortality.single.run , lifespan=tb$rls, 
                method="L-BFGS-B", lower=c(1E-10,0.05,1E-1), upper=c(10,100,50) );  
  
  output[i, c("R", "t0", "n")] = fitBinom$par
    
  tb2 = calculate.mortality.rate(tb$rls)
  plot( tb2$s ~ tb2$t, main=filename)
  plot( tb2$mortality.rate ~ tb2$t , main=filename)
  plot( tb2$mortality.rate ~ tb2$t, log='y', main='linear for Gompertz' )
  plot( tb2$mortality.rate ~ tb2$t, log='yx', main='linear for Weibull' )
}#i file loop

output = output[order(output$cycloheximide), ]

#plot( output$R ~ output$cycloheximide, type='l' )
plot( output$t0 ~ output$cycloheximide,type='p', xlab="cycloheximide", ylab="t0" )
points( output$cycloheximide, output$t0, pch=19)
m= lm(output$t0 ~ output$cycloheximide ); abline(m, col='red'); summary(m)

plot( output$n ~ output$cycloheximide,type='l',  xlab='cycloheximide', ylab="n" )
points( output$cycloheximide, output$n, pch=19)
#m= lm(output$n ~ output$cycloheximide ); abline(m, col='red')

plot( output$mean ~ output$cycloheximide, type='p', xlab="cycloheximide",ylab="Mean RLS")
points( output$cycloheximide, output$mean, pch=19)
m= lm(output$mean ~ output$cycloheximide ); abline(m, col='red')
#plot( output$CV ~ output$cycloheximide, type='l', xlab="cycloheximide", ylab="CV of RLS" )
#points(output$cycloheximide, output$CV, pch=19)

#plot( output$stddev ~ output$cycloheximide, type='l')

plot( output$Rflex ~ output$cycloheximide, type='p', xlab="cycloheximide",ylab="R")
points( output$cycloheximide, output$Rflex, pch=19)
plot( output$Gflex ~ output$cycloheximide, type='p', xlab="cycloheximide",ylab="G")
points( output$cycloheximide, output$Gflex, pch=19)

dev.off()


q('no')
####END of 20140712 #########


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

