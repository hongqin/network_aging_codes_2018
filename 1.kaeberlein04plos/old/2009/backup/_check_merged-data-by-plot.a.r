#092304Thu
 
# jpeg(  "092104_summary_survival_curves.jpg", width=480, height=480,
#        pointsize=12, quality=600, bg="white" );
postscript ( "092304.mine.vs.kaeberlein.ps", width=8, height=8 )

 library(survival);
 rls <- read.table("rls.tab", sep="\t", header=T);
 strains <- names(rls);
 
  num.col <- length( as.list( rls[1,]) )
 num.row <- length( as.list( rls[,1]) )
 cols <- rainbow(num.col);
 
 ylim <- c(0,1);
 xlim <- c(0,100);
 
 tmp.fit <- survfit ( Surv(rls[,1]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 	);
 
 for( i in 2 : num.col) {
  tmp.fit <- survfit( Surv(rls[,i]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i] );
 }
  
 leg.txt <- names(rls)
 legend( 60, 0.9, leg.txt, col=cols, lty=1,lwd=2 );

 
##kaeberlein fig 1a 
 start = 10; end = 14;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 	);
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##fig 1a end
  
  
 
 
##kaeberlein fig 1b
 start = 15; end = 18;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 	);
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
## end 


 
##kaeberlein fig 1c
 start = 19; end = 22;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 



 
##kaeberlein fig 2a
 start = 23; end = 26;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 


 
##kaeberlein fig 2b
 start = 27; end = 29;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 



 
##kaeberlein fig 2c
 start = 30; end = 32;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 



 
##kaeberlein fig 3
 start = 36; end = 43;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 

 ks.test( rls[,43], rls[,41]); #p-value = 0.9 similar in t.test 
 #Yes. Now, that's why Kaeberlein use bar plot.
 
 
 

 
##kaeberlein fig 4a
 start = 44; end = 46;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 


 

 
##kaeberlein fig 4b
 start = 47; end = 50;  
 cols <- rainbow( end - start + 1);
 max.rls = max( rls[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(rls[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 );
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(rls[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(rls[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
##end 



 
 dev.off();

 quit("yes");

