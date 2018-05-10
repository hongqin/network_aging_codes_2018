#042305Sat

library(nlme); library(survival);

rls = read.table( "rls.tab", header=T );
expts = names(rls);
expts;

####### 
i = 1;
tmp.fit <- survfit( Surv(rls[,i]) )
expts[i];
tmp.fit

####### 
i = 2;
tmp.fit <- survfit( Surv(rls[,i]) )
expts[i];
tmp.fit

####### 
i = 3;
tmp.fit <- survfit( Surv(rls[,i]) )
expts[i];
tmp.fit
  
####### 
i = 4;
tmp.fit <- survfit( Surv(rls[,i]) )
expts[i];
tmp.fit


quit( "yes" );

