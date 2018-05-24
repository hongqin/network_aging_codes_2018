#test binomial aging model on SIR2

rm(list=ls());
setwd("~/projects/0.network.aging.prj/0a.network.fitting")
source("lifespan.20140711.r")
require(flexsurv)

#files = list.files(, path='DR', pattern="4742")
files = c( 
  "BY4742_MATalpha_temp30_media_YPD_n19930.csv",     #WT
  "sir2_MATalpha_temp30_media_YPD_n1497.csv",        #sir2
  "SIR2OE_MATalpha_temp30_media_0.05% D_n240.csv",   #SIR2
  "SIR2OE_MATa_temp30_media_YPD_n341.csv"        
);

i=1; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s1 = calculate.s(tb$rls); s1=rbind(c(1,0), s1)
i=2; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s2 = calculate.s(tb$rls); s2=rbind(c(1,0),s2);
i=3; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s3 = calculate.s(tb$rls); s3=rbind(c(1,0),s3);
i=4; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s4 = calculate.s(tb$rls); s4=rbind(c(1,0),s4);

plot( s1$s ~ s1$t, xlim=c(1,100), type='l', col='black', )
lines(s2$s ~ s2$t, col='red')
lines(s3$s ~ s3$t, col='blue')
lines(s4$s ~ s4$t, col='orange')

output = data.frame(files); 
output$mean = NA; output$stddev = NA; output$CV=NA; 
output$R=NA; output$t0=NA; output$n=NA;
output$Rflex = NA; output$Gflex = NA; 

pdf("_sir2-binomial-20140801.pdf", width=10, height=8)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))


#i=1; 
for( i in 1:3) {
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
                method="L-BFGS-B", lower=c(1E-10,5,1), upper=c(10,100,50) );  
  
  output[i, c("R", "t0", "n")] = fitBinom$par
  
    #  fitGomp2 = optim ( c(0.005, 0.08),  llh.G.single.run , lifespan=tb$rls,
   #                    method="L-BFGS-B", lower=c(1E-10,1E-2), upper=c(10,10) );  
  
  tb2 = calculate.mortality.rate(tb$rls)
  plot( tb2$s ~ tb2$t, main=filename)
  plot( tb2$mortality.rate ~ tb2$t , main=filename)
  plot( tb2$mortality.rate ~ tb2$t, log='y', main='linear for Gompertz' )
  plot( tb2$mortality.rate ~ tb2$t, log='yx', main='linear for Weibull' )
}#i file loop

dev.off()

#"BY4742_MATalpha_temp30_media_YPD_n19930.csv",     #WT
#"sir2_MATalpha_temp30_media_YPD_n1497.csv",        #sir2
#"SIR2OE_MATalpha_temp30_media_0.05% D_n240.csv",   #SIR2

i=1; filename = paste('DR/', files[i], sep='');  WT = read.csv(filename)
i=2; filename = paste('DR/', files[i], sep='');  sir2 = read.csv(filename)
i=3; filename = paste('DR/', files[i], sep='');  SIR2 = read.csv(filename)

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

########### LRT on sir2 and WT
llh.H0  = optim( c(0.01,0.1,0.01, 0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=sir2, lifespan2=WT );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=sir2, lifespan2=WT );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=sir2, lifespan2=WT );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=sir2, lifespan2=WT );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);LH
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );
deltaLH =  - LH + llh.H1i$value; 
1 - pchisq( 2*deltaLH, df =1 );


########### LRT on SIR2 and WT
llh.H0  = optim( c(0.01,0.1,0.01, 0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=SIR2, lifespan2=WT );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=SIR2, lifespan2=WT );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=SIR2, lifespan2=WT );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=SIR2, lifespan2=WT );
cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);LH
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );
deltaLH =  - LH + llh.H1i$value; 
1 - pchisq( 2*deltaLH, df =1 );


#######LRT on sir2 and SIR2
llh.H0  = optim( c(0.01,0.1,0.01,0.1), Hx.llh.gompertz.single.run, model=H0,  lifespan1=sir2, lifespan2=SIR2 );
llh.H1i = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1i, lifespan1=sir2, lifespan2=SIR2 );
llh.H1g = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H1g, lifespan1=sir2, lifespan2=SIR2 );
llh.H2  = optim( c(0.01,0.2,0.01,0.1), Hx.llh.gompertz.single.run, model=H2,  lifespan1=sir2, lifespan2=SIR2 );

cbind(llh.H0$par, llh.H1i$par, llh.H1g$par, llh.H2$par);
LH = c(llh.H0$value, llh.H1i$value, llh.H1g$value, llh.H2$value);LH
deltaLH =  - LH + llh.H0$value; 
1 - pchisq( 2*deltaLH, df =1 );
deltaLH =  - LH + llh.H1i$value; 
1 - pchisq( 2*deltaLH, df =1 );

###################



#q('no')
####END of 20140801 #########


