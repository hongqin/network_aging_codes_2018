## the main section 
 tb <- read.table( "wm.out.tab", header=T);

 x = cbind( mean(tb$a), sd(tb$a), mean(tb$b), sd(tb$b));

 boot <- read.table( "boot.tab", header=F);
 rls = mean( boot );

 x = cbind( mean(rls), sd(rls), x );
 
 x = round( x, 7);

 out = data.frame( x );
 names( out ) = cbind("rls", "rls.sd", "a", "a.sd", "b", "b.sd");
 out;

 write.table( out, "summary.wm.out.tab", row.names=F, col.names=T, quote=F, sep="\t");

q("yes")
