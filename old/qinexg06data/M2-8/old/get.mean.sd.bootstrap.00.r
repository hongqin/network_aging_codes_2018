#040405Mon

tb = read.table( "bootstrap.rls.tab", header=F );
x = mean(tb)

avg.rls = mean( x );
sd.rls = sd( x );

### mean rls
avg.rls
### sd rls
sd.rls

quit( "yes" );

