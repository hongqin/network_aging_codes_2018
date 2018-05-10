rm(list=ls());

filesNlabels = c( 
"010305.BY4743.rls.tab", "by4743a", #1
"011705.BY4743.rls.tab", "by4743b", #3
"122004.BY4743.rls.tab", "by4743c", #5
"050905.BY4742.rls.tab", "by4742a", #7
"053104.BY4742.rls.tab", "by4742b", #9
"061004.BY4742.rls.tab", "by4742c", #11
"053104.BY4741.rls.tab", "by4741a", #13
"061004.BY4741.rls.tab", "by4741b", #15
"042805.BY4716.rls.tab", "by4716",  #17
"050905.sir2D.4742.rls.tab", "sir2D4742a",  #19
"051805.sir2D.4742.rls.tab", "sir2D4742b",  #21
"052604.sir2D.4741a.rls.tab", "sir2D4741",  #23
"050905.sir2DSIR2.4742.rls.tab","sir2DSIR24742a", #25
"051805.sir2DSIR2.4742.rls.tab","sir2DSIR24742b", #27
"062404.sir2D.4741a.CR.rls.tab","sir2DCR4741a",   #29
"062904.by4741.CR.rls.tab","by4741CR");           #31

files  = filesNlabels[seq(1,31, by=2) ];
labels = filesNlabels[seq(2,32, by=2) ];

tb = read.table( files[1] );
len = 50;  #lenth of column
tb[ len,1] =  NA;  #add an extral elements
for( i in 2:length(files) ) {
  tmp = read.table( files[i] );
  del = len - length( tmp[,1] );
  if ( del > 0 ) {
   tb = cbind( tb, c(tmp[,1], rep(NA, del) ) );
  } else {
  tb = cbind( tb, tmp );
  }
} 
colnames(tb) = labels;

##### log likelihood function, simple gompertz mortality model
   #s = exp( (I/G) *(1 - exp(G* my.data)) )  ;
   #m = I * exp( G * my.data );   
 llh.gompertz.single.run <- function( IG, lifespan ) {
   I = IG[1]; G = IG[2]; 
   my.data = lifespan[!is.na(lifespan)];
   log_s = (I/G) *(1 - exp(G* my.data))
   log_m = log(I) +  G * my.data ; 
   my.lh = sum(log_s)  + sum(log_m);
   print (IG ); #trace the convergence
   ret = - my.lh # because optim seems to maximize 
 }

#####

x1 = llh.gompertz( 0.004, 0.08);
x2 = llh.gompertz( 0.04, 0.08);

ret1 = optim ( c(0.005, 0.08),  llh.gompertz.single.run, lifespan=tb[,1] );
ret2 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=tb[,1] );  
# nice, converge to the same values, same as Qin06EXG publication!

fit = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=tb[,3] );  

### partition and fitting by strains
by4743 = c( tb[,1], tb[,2], tb[,3]);
fit4743 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4743 );  

by4742 = c( tb[,4], tb[,5], tb[,6]);
fit4742 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4742 );  

#by4741 = c( tb[,7], tb[,8], tb[,9]);
by4741 = c( tb[,7], tb[,8] );
fit4741 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4741 );  

by4716 = c( tb[,10]);
fit4716 = optim ( c(0.01, 0.2),  llh.gompertz.single.run, lifespan=by4716 );  


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

