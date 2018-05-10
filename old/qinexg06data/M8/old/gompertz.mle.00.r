#040405Mon

library(nlme); library(survival);

###Gompertzian survival model
# surival chance
gomp.s <- function( I, G, t ) {
#   tmp1 = 1 - exp( G * t );
   ret <- exp( (I/G) *(1 - exp(G*t)) )  ;
}

gomp.3p.s <- function( I,G,M, t) {
# tmp1 = -M*t
 #tmp2 = (I/G)
 exp( -M * t + (I/G)*(1- exp(G*t)))
}

#mortality rate
gompertiz.mortality.rate <- function( I, G, t) { ret <- I * exp(G *t) } 
gomp.3p.m <- function( I,G,M,t) { I * exp(G*t) + M }

#log desity function of death
log.density.death <- function( age, I, G ) {
   s = gompertz.viability ( I, G, age);
   m = gompertiz.mortality.rate (I, G, t) 
   ret <-  log(s) + log( m )
}

#likelihood function by density
lh.gomp <- function( data, I, G ) {
   my.data = data[!is.na(data)];
   s = exp( (I/G) *(1 - exp(G* my.data)) )  ;    # print( s );
   m = I * exp( G * my.data );    # print( m);
   my.lh = log(s) + log(m);
   ret <- sum( my.lh);
}

#likelihood function by survival function
lh.gomp.s <- function( data, I, G ) {
   my.data = data[!is.na(data)];
   #s = exp( (I/G) *(1 - exp(G* my.data)) )  ;     print( s );
   s1 = exp( (I/G) *(1 - exp(G* (my.data-0.5))) )  ;     # print( s1 );
   s2 = exp( (I/G) *(1 - exp(G* (my.data+0.5))) )  ;     # print( s2 );
   my.lh = log(s1 - s2)
   ret <- sum( my.lh);
}

###end of Gompertzain model###

tb = read.table( "bootstrap.rls.tab", header=F );
col.num.in = length( tb[1,] );
tb.mle.out.d = data.frame( matrix(ncol=5, nrow=col.num.in ) );
names( tb.mle.out.d ) = c("I.mle", "I.0", "G.mle", , "G.0", "lh.max" );
tb.mle.out.s = tb.mle.out.d;

for ( k in 1:col.num.in ) {
 rls = tb[ , k ];
 # calculate the initial values
 tmp.fit = survfit( Surv(rls) );
 s = c( 1, tmp.fit$surv );
 t = c( 0, tmp.fit$time );
 fm <- gnls( s ~ exp( (I/G)* ( 1 - exp(G * t))),start = list( I = 0.001, G = 0.2));
 I.0 <- fm$coefficients[1];
 G.0 <- fm$coefficients[2];  

 lh.gomp.0 = lh.gomp( rls, I.0, G.0);
 lh.gomp.0;
 lh.gomp.0s = lh.gomp.s( rls, I.0, G.0);
 lh.gomp.0s;


 #define the search space 
 I.range = c( 0, min(10 * I.0, 0.1) );
 G.range = c( 0, min( 3* G.0, 0.5 ) );
 I.itr = 1000;
 G.itr = 200;
 I.step = ( I.range[2] - I.range[1] ) / I.itr;
 G.step = ( G.range[2] - G.range[1] ) / G.itr; 

 Is = seq( I.step, I.range[2], length= I.itr );
 Gs = seq( G.step, G.range[2], length= G.itr );
 lh.mat = matrix( data=NA, nrow=I.itr, ncol=G.itr );

 lh.max = -Inf;
 I.max = NA;
 G.max = NA;

 for ( i in 1:I.itr  ) {
   for ( j in 1:G.itr  ) {
     lh.mat[i,j] = lh.gomp( rls, Is[i], Gs[j] );
     if ( lh.mat[i,j] >= lh.max ) {
	 lh.max = lh.mat[i, j];
	 I.max = i;
	 G.max = j;
     } 
   }
 }

 max( lh.mat );
 lh.max; I.max; G.max;
 tb.mle.out.d [k, ] = c(Is[I.max], I.0, Gs[G.max], G.0, lh.max ); ####

 print ( k ); 

 #Now try mle by survival function

 lh.mat.s = matrix( data=NA, nrow=I.itr, ncol=G.itr );
 lh.max.s = -Inf;
 I.max.s = NA;
 G.max.s = NA;

 for ( i in 1:I.itr  ) {
  for ( j in 1:G.itr  ) {
    lh.mat.s[i,j] = lh.gomp.s( rls, Is[i], Gs[j] );
    if ( lh.mat.s[i,j] >= lh.max.s ) {
	lh.max.s = lh.mat.s[i, j];
	I.max.s = i;
	G.max.s = j;
    } 
  }
  }

 max( lh.mat.s );
 c(lh.max.s, I.max.s, G.max.s, Is[I.max.s], I.0, Gs[G.max.s], G.0); 
 c(lh.max, I.max, G.max, Is[I.max], I.0, Gs[G.max], G.0); 

 tb.mle.out.s [k, ] = c(Is[I.max.s], I.0, Gs[G.max.s], G.0, lh.max.s ); ####

 print ( k ); 
} # end of k loop

write.table( tb.mle.out.d, "den-mle.out.bootstrap.rls.tab", sep="\t", row.names=F,col.name=T );
write.table( tb.mle.out.s, "sur-mle.out.bootstrap.rls.tab", sep="\t", row.names=F,col.name=T );

### mle by density results
summary( tb.mle.out.d );
sd( tb.mle.out.d );

### mle by survival curve  results
summary( tb.mle.out.s );
sd( tb.mle.out.s );

quit( "yes" );

