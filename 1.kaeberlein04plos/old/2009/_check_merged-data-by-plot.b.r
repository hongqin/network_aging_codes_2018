#092304Thu
 
# jpeg(  "092104_summary_survival_curves.jpg", width=480, height=480,
#        pointsize=12, quality=600, bg="white" );
# postscript ( "092304.mine.vs.kaeberlein.ps", width=8, height=8 )

 library(survival);
 rls <- read.table("rls.tab", sep="\t", header=T);
 strains <- names(rls);
 ylim = c( 0,1);
 
  num.col <- length( as.list( rls[1,]) )
 num.row <- length( as.list( rls[,1]) )
 
### how many BY4742 data set in this paper?
 tmp = strains[ grep ("4742", strains)  ]
 tmp2 = tmp[ grep( "fig", tmp) ];
 tmp2[13:14] = NA;
 by4742.used.in.Kaeberlein04Plos = tmp2[! is.na(tmp2) ]
 
 tb.mkby4742 = rls[, by4742.used.in.Kaeberlein04Plos ];
  
 postscript("BY4742.in.Kaeberlein04Plos.ps");

 pairs( tb.mkby4742 ); 
 #this shows that only fig 2c and fig 2d use the same BY4742 data
 # there are many linar correlation in the initial part, but leveled off in the end.
 
 start = 1; end = length( names(tb.mkby4742)) ;  
 cols <- rainbow( end - start + 1);
 max.rls = max( tb.mkby4742[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(tb.mkby4742[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 	);
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(tb.mkby4742[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(tb.mkby4742[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
 dev.off();
### how many BY4742 data set in this paper?  END
 
 
 
### how many sir2 data set in this paper?
 tmp = strains[ grep ("sir2", strains)  ]
 tmp2 = tmp[ grep( "fig", tmp) ];
 #tmp2[c(14)] = NA;
 sir2.used.in.Kaeberlein04Plos = tmp2[! is.na(tmp2) ]
 
 tb.mksir2 = rls[, sir2.used.in.Kaeberlein04Plos ];
  
 postscript("sir2.in.Kaeberlein04Plos.ps");

 pairs( tb.mksir2 ); 
 #this shows that only fig 2c and fig 2d use the same BY4742 data
 # there are many linar correlation in the initial part, but leveled off in the end.
 
 start = 1; end = length( names(tb.mksir2)) ;  
 cols <- rainbow( end - start + 1);
 max.rls = max( tb.mksir2[,start:end], na.rm=T )
 xlim = c( 0,  max.rls ); 
 tmp.fit <- survfit ( Surv(tb.mksir2[, start]) );
 y <- c( 1, tmp.fit$surv );
 x <- c( 0, tmp.fit$time );
 plot 	( x, y, col=cols[1], ylim=ylim, xlim=xlim, type='l', 
	xlab="generations", ylab="Viability",
        main="H Qin Sep 23, 2004"
 	);
 
 for( i in 1:(end - start) ) {
  tmp.fit <- survfit( Surv(tb.mksir2[, (i + start) ]) )
  y <- c( 1, tmp.fit$surv );
  x <- c( 0, tmp.fit$time );
  lines ( x, y, col=cols[i+1] );
 }
  
 leg.txt <- names(tb.mksir2[,start:end])
 legend( (max.rls * 0.70 ), 0.9, leg.txt, col=cols, lty=1,lwd=2 );
 dev.off();
### how many sir2 data set in this paper?  END

 quit("no");

