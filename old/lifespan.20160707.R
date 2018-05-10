# source("/Users/hongqin/lib/R/lifespan.r")

calculate.mortality.rate = function( lifespan ) {
  tb = calculate.s(lifespan)
  tb$ds=NA; tb$dt=NA
  
  tb$dt[1] = NA; 
  tb$ds[1] = NA; 
  #tb$mortality.rate[1] = tb$ds[1] / tb$dt[1]
  tb$mortality.rate[1] = NA
  
  for( j in 2:length(tb[,1])) {
    tb$ds[j] =  tb$s[j-1] - tb$s[j] 
    tb$dt[j] = -tb$t[j-1] + tb$t[j]
    tb$mortality.rate[j] = tb$ds[j] / ( tb$s[j] * tb$dt[j])
  }
  return( tb )
}

##### log likelihood function, 3-parameter binomial mortality rate model
# http://hongqinlab.blogspot.com/2013/12/binomial-mortailty-model.html
# m = R ( 1 + t/t0)^(n-1)
# s = exp( (R t0/n)*(1 - (1+t/t0)^n ) )  
# G = (n-1)/t0

llh.binomialMortality.single.run <- function( Rt0n, lifespan, debug=0 ) {
  #to do: check make sure Rt0n has the right format. 
  I = Rt0n[1]; t0 = Rt0n[2]; n=Rt0n[3];
  my.data = lifespan[!is.na(lifespan)];
  log_s = (I * t0 /n )*(1 - (1 + my.data/t0)^n);
  log_m = log(I) +  (n-1) * log(1 + my.data/t0 ); 
  my.lh = sum(log_s)  + sum(log_m);
  if(debug>0) { print (paste("llh.binomialMortality.single.run:", t(Rt0n) ))}  #trace the convergence
  ret = - my.lh # because optim seems to maximize 
}

#random number from binomial-aging model
rbinomial.aging <- function ( Rt0n, n){
  x.uniform = runif(n)
  inverse.binomialaging.CDF = function(R,t0,n,y) {t0*((1- n*log(y)/(R*t0))^(1/n) -1) }
  x.binomialaging = inverse.binomialaging.CDF(Rt0n[1],Rt0n[2], Rt0n[3], x.uniform)
  return(x.binomialaging)
}

##########################################
#inverse of gompertz CDF
# see http://hongqinlab.blogspot.com/2013/06/median-lifespan-of-2-parameter-gompertz.html
#inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }

#see generate random number with a given distribution
# http://hongqinlab.blogspot.com/2013/12/generate-gompertz-random-numbers.html

#generate Gompertz random numbers
rgompertz = function(R,G, n){
  x.uniform = runif(n)
  inverse.gomp.CDF = function(R,G,y) {  (1/G)*log(1 - (G/R)*log(1-y)  ) }
  x.gompertz = inverse.gomp.CDF(R,G, x.uniform) #20160711 corret an error here. 
  return(x.gompertz)
}
rgompertz(0.001,0.2,100)

##########################################
####
 #calculate s
 calculate.s = function( lifespan ){
 	 my.data = sort( lifespan[!is.na(lifespan)] );
   tmptb = table( my.data )
    for( i in 2:length(tmptb)) {
    	tmptb[i] = tmptb[i-1] + tmptb[i]    		} 	 
 	  tot = length(my.data)
 	 tmptb = tmptb / tot; 
 	 s = 1 - tmptb
 	 #list( s=s, t=unique(my.data));
 	 ret = data.frame( cbind(s, unique(my.data)));
 	 names(ret) = c("s", "t");
 	 ret;
 	}

##########################################
##### log likelihood function, simple gompertz mortality model
   #s = exp( (I/G) *(1 - exp(G* my.data)) )  ;
   #m = I * exp( G * my.data ) ;   
 llh.G.single.run <- function( IG, lifespan ) {
   I = IG[1]; G = IG[2];  
   my.data = lifespan[!is.na(lifespan)];
   log_s = (I/G) *(1 - exp(G* my.data)) 
   #log_m = log( I * exp(G*my.data) + M); 
   log_m = log( I) + G*my.data; 
   my.lh = sum(log_s)  + sum(log_m);
   print (c(IG, my.lh) ); #trace the convergence
   ret = - my.lh # because optim seems to minimize
 }

##### log likelihood function, Gompertz-Makeham mortality model 
   #s = exp( (I/G) *(1 - exp(G* my.data)) - M*t )  ;
   #m = I * exp( G * my.data ) + M;   
 llh.GM.single.run <- function( IGM, lifespan ) {
   I = IGM[1]; G = IGM[2]; M = IGM[3]; 
   my.data = lifespan[!is.na(lifespan)];
   log_s = (I/G) *(1 - exp(G* my.data)) - M*my.data
   log_m = log( I * exp(G*my.data) + M); 
   #log_m = log( I) + G*my.data; 
   my.lh = sum(log_s)  + sum(log_m);
   print (c(IGM, my.lh) ); #trace the convergence
   ret = - my.lh # because optim seems to minimize
 }

##########################################
   #s = exp( (I/G) *(1 - exp(G* my.data)) - M*t )  ;
 GM.s = function( IGM, t ) {
 	 I = IGM[1]; G = IGM[2]; M = IGM[3]; 
   log_s = (I/G) *(1 - exp(G* t)) - M*t
   ret  = exp( log_s )
 	}

#####
 G.s = function( IGM, t ) {
   I = IGM[1]; G = IGM[2]; M = 0; 
   log_s = (I/G) *(1 - exp(G* t)) - M*t
   ret  = exp( log_s )
 }
 
 #####
 
