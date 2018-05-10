#092304Thu

## the main section 
 tb <- read.table( "wm.out.tab", header=T);

 x = cbind( mean(tb$a), sd(tb$a), mean(tb$b), sd(tb$b));
 x = round( x, 6);

 out = data.frame( x );
 names( out ) = cbind("a", "a.sd", "b", "b.sd");
 out;

 write.table( out, "summary.wm.out.tab", row.names=F, col.names=T, quote=F, sep="\t");

q("yes")
