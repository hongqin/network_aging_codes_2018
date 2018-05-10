rm (list=ls())
setwd("~/projects/network.aging.prj/lnR-G/qinexg06data")

list.files()
#tb = read.table("061405.rls.for.publication.tab", header=T)
tb = read.table("052305.strain.a.sd.b.sd.tab", header=T)
tb$strain = as.character(tb$strain)
row.names(tb) = tb$strain
tb.ori = tb

my2Nstrains = c("101S",  "BY4743", "M1-2", "M13", "M14", "M2-8", "M22", "M32", "M34", "M5", "M8", "RM112N", 
  "RXB", "SGU57", "SK1",   "W303", "YPS128", "YPS163")

myNatstrains = c("101S",   "M1-2", "M13", "M14", "M2-8", "M22", "M32", "M34", "M5", "M8", "RM112N", "SGU57", "YPS128", "YPS163")

tb = tb.ori[myNatstrains, ] #p=0.065
summary(lm( log(tb$a)~ tb$b ),  )

tb = tb.ori[my2Nstrains, ] #p=0.042
summary(lm( log(tb$a)~ tb$b ),  )
m = lm(log(tb$a) ~ tb$b)

tiff("lnR-G-diploid-20130617.tif",width=480,heigh=480)
par(font=2)
plot( log(tb$a) ~ tb$b, pch=19, col='blue', xlim=c(0.05, 0.22),ylim=c(-8,-3.5),xlab="G",ylab="lnR") 
text( tb$b*1.02, log(tb$a*1.1), as.character(tb$strain ) )
abline(m, col='red')
text(0.15, -4.1, "R^2=0.23, p=0.042")
dev.off()


