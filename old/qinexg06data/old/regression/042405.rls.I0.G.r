library(ade4)

tb.all = read.table( "merged.tab", header=T);
expt.name = as.character(tb.all$expt);

strain = expt.name;
for( i in 1:length(expt.name)) {
   strain[i] = unlist( strsplit( expt.name[i], "\\.") ) [1];
}
strain;

plot( avg~a, data=tb.all, xlab="I0", ylab="Average RLS",
 main= "Average lifespan ~ I0, merged, WinModest, 042405 " );
 
scatterutil.eti( tb.all$a, tb.all$avg, strain, 0.5 );
#arrows( tb.all$a_low, tb.all$avg, tb.all$a_up, tb.all$avg, col="green", angle=90, lty=8, length=0.1,code=3  );
# 95% confidence level is too stringent

plot( avg~b, data=tb.all, xlab="G", ylab="Average RLS",
 main= "Average lifespan ~ G, merged, WinModest, 042405 " );
scatterutil.eti( tb.all$b, tb.all$avg, strain, 0.5 );

###############################################
# regression without BY4743
tb.all[2,]=NA;
rls.avg.nat = mean( tb.all$avg, na.rm=T );

m.a = lm( avg ~ a, data=tb.all );
summary(m.a);
a.avg = mean( tb.all$a, na.rm=T );

m.b = lm( avg ~ b, data=tb.all );
summary(m.b);
b.avg = mean( tb.all$b, na.rm=T  );

##############################################
##### theoretic predictions 
avg.predict <- function ( a, b ) {   log( log(2)*b/a + 1 )/b  }

x.a = seq(0.00001, max(tb.all$a, na.rm=T), length=500);
y.rls.a = avg.predict( x.a, b.avg);

x.b = seq( 0.05, max(tb.all$b, na.rm=T), length=500);
y.rls.b = avg.predict( a.avg, x.b );


tb.all = read.table( "merged.tab", header=T);
expt.name = as.character(tb.all$expt);

strain = expt.name;
for( i in 1:length(expt.name)) {
   strain[i] = unlist( strsplit( expt.name[i], "\\.") ) [1];
}
strain;

####################################3
### redo plot with prediction and regression
plot( avg~a, data=tb.all, xlab="I0", ylab="Average RLS",
 main= "Average lifespan ~ I0, merged, WinModest, 042405 " ); 
scatterutil.eti( tb.all$a, tb.all$avg, strain, 0.5 );
abline( m.a, col="red", lty=2 );
lines( x.a, y.rls.a, col="blue", lty=2);

#arrows( tb.all$a_low, tb.all$avg, tb.all$a_up, tb.all$avg, col="green", angle=90, lty=8, length=0.1,code=3  );
# 95% confidence level is too stringent

plot( avg~b, data=tb.all, xlab="G", ylab="Average RLS",
 main= "Average lifespan ~ G, merged, WinModest, 042405 " );
scatterutil.eti( tb.all$b, tb.all$avg, strain, 0.5 );
abline( m.b, col="red", lty=2 );
lines( x.b, y.rls.b, col="blue", lty=2);
