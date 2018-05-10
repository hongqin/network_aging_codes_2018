calculate.mortality.rate = function( lifespan ) {
  tb = calculate.s(lifespan)
  tb$ds=NA; tb$dt=NA
  
  tb$dt[1] = tb$s[1]
  tb$ds[1] = 1 - tb$s[1]
  tb$mortality.rate[1] = tb$ds[1] / tb$dt[1]
  
  for( j in 2:length(tb[,1])) {
    tb$ds[j] =  tb$s[j-1] - tb$s[j] 
    tb$dt[j] = -tb$t[j-1] + tb$t[j]
    tb$mortality.rate[j] = tb$ds[j] / ( tb$s[j] * tb$dt[j])
  }
  return( tb )
}