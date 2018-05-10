#test binomial aging model on glucose dose in BY4742

rm(list=ls());
setwd("~/github/0.network.aging.prj/0a.network.fitting")
source("lifespan.20140711.r")
require(flexsurv)

#files = list.files(, path='DR', pattern="BY4742")
files = c(
  "BY4742_MATalpha_temp30_media_0.005% D_n299.csv",                   
  "BY4742_MATalpha_temp30_media_0.05% D_n4690.csv",                   
  "BY4742_MATalpha_temp30_media_0.1% D_n338.csv",                     
  "BY4742_MATalpha_temp30_media_0.5% D_n579.csv",                     
  "BY4742_MATalpha_temp30_media_YPD_n19930.csv"
);

i=1; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s1 = calculate.s(tb$rls); s1=rbind(c(1,0), s1)
i=2; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s2 = calculate.s(tb$rls); s2=rbind(c(1,0),s2);
i=3; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s3 = calculate.s(tb$rls); s3=rbind(c(1,0),s3);
i=4; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s4 = calculate.s(tb$rls); s4=rbind(c(1,0), s4);
i=5; filename = paste('DR/', files[i], sep='');  tb = read.csv(filename)
s5 = calculate.s(tb$rls); s5=rbind(c(1,0), s5)

plot( s1$s ~ s1$t, xlim=c(10,40), type='l', col='red', )
lines(s2$s ~ s2$t, col='orange')
lines(s3$s ~ s3$t, col='blue')
lines(s4$s ~ s4$t, col='green')
lines(s5$s ~ s5$t, col='black')


output = data.frame(files); 
output$mean = NA; output$stddev = NA; output$CV=NA; 
output$R=NA; output$t0=NA; output$n=NA;
output$Rflex = NA; output$Gflex = NA; 
output$glucose = c(0.005, 0.05, 0.1, 0.5, 2.0)


#pdf("BY4742GlucoseDose-binomial-20140713.pdf", width=10, height=8)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))


#i=1; 
for( i in 1:5) {
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
  
    #  fitGomp2 = optim ( c(0.005, 0.08),  llh.G.single.run , lifespan=tb$rls,
   #                    method="L-BFGS-B", lower=c(1E-10,1E-2), upper=c(10,10) );  
  
  tb2 = calculate.mortality.rate(tb$rls)
  plot( tb2$s ~ tb2$t, main=filename)
  plot( tb2$mortality.rate ~ tb2$t , main=filename)
  plot( tb2$mortality.rate ~ tb2$t, log='y', main='linear for Gompertz' )
  plot( tb2$mortality.rate ~ tb2$t, log='yx', main='linear for Weibull' )
}#i file loop


plot( output$R ~ output$glucose, type='l' )
plot( output$t0 ~ output$glucose,type='l', xlab="Glucose", ylab="t0" )
points( output$glucose, output$t0, pch=19)

plot( output$n ~ output$glucose,type='l', ylim=c(4.4,5.5), xlab='Glucose', ylab="n" )
points( output$glucose, output$n, pch=19)

plot( output$mean ~ output$glucose, type='l', xlab="Glucose",ylab="Mean RLS")
points( output$glucose, output$mean, pch=19)

plot( output$CV ~ output$glucose, type='l', xlab="Glucose", ylab="CV of RLS", ylim=c(0.35, 0.39) )
points(output$glucose, output$CV, pch=19)

plot( output$stddev ~ output$glucose, type='l')

plot( output$Rflex ~ output$glucose, type='l', xlab="Glucose",ylab="R")
points( output$glucose, output$Rflex, pch=19)
plot( output$Gflex ~ output$glucose, type='l', xlab="Glucose",ylab="G")
points( output$glucose, output$Gflex, pch=19)

#dev.off()


#q('no')
####END of 20140712 #########
