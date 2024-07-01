rm(list=ls())
library(spatstat)

#
# Location of male, female, and juvenile tupelo trees in 3 plots in the Savannah River tupelo
# Courtesy of Philip Dixon
#
tupelo<-read.table("tupelo.txt",header=TRUE)
tupelo<-tupelo[tupelo$sex != "j",]

# number of plot
n.plot<-3
tupelo.ppp<-ppp(tupelo[tupelo$plot == n.plot ,5],tupelo[tupelo$plot == n.plot,6], marks=tupelo$sex[tupelo$plot == n.plot], window =   owin(xrange=c(-1,53), yrange=c(0,53)))

plot(tupelo.ppp$x,tupelo.ppp$y, xlab='x',ylab='y',type='n')
points(tupelo.ppp$x[tupelo.ppp$marks =='m'],tupelo.ppp$y[tupelo.ppp$marks =='m'], pch=20,cex=1.2)
points(tupelo.ppp$x[tupelo.ppp$marks =='f'],tupelo.ppp$y[tupelo.ppp$marks =='f'], pch='+',cex=1.2)




# Random labelling approach
# we select the male

X.m <- unmark(split(tupelo.ppp)$m)
K.m<-Kest(X.m)
int.m<-summary(X.m)$intensity
# we select the female

X.f <- unmark(split(tupelo.ppp)$f)
K.f<-Kest(X.f)
int.f<-summary(X.f)$intensity
K.mf<-Kcross(tupelo.ppp,"m","f")
K.fm<-Kcross(tupelo.ppp,"f","m")
K<-(int.m*K.mf$iso+int.f*K.fm$iso)/(int.m+int.f)
diff<-K.m$iso-K.f$iso
diff2<-K.f$iso-K
r<-K.mf$r
m<-40
n<-length(diff)
sdiff<-matrix(0,m,n)
sdiff2<-sdiff

set.seed(19)
for (i in 1:m) {
      smarks<-sample(tupelo.ppp$marks)
      X<-tupelo.ppp
      X$marks<-smarks	
      X.m<-unmark(split(X)$m)
      K.m<-Kest(X.m)
      X.f<-unmark(split(X)$f)
      K.f<-Kest(X.f)
      K.mf<-Kcross(X,"m","f")
      K.fm<-Kcross(X,"f","m")
      sdiff[i,]<-K.m$iso-K.f$iso
      K<-(int.m*K.mf$iso+int.f*K.fm$iso)/(int.m+int.f)
      sdiff2[i,]<-K.f$iso-K

}
sll<-apply(sdiff,2,min)
sul<-apply(sdiff,2,max)


plot(K.f$r,diff,type='l',lty=1,xlab="h", ylab= expression(hat(D)(h)), ylim=c(min(sll),max(sul)))
lines(K.f$r,sll,type='l',lty=2)
lines(K.f$r,sul,type='l',lty=2)



sll<-apply(sdiff2,2,min)
sul<-apply(sdiff2,2,max)

plot(r,diff2,type='l',lty=1,xlab="h", ylab= expression(hat(D)[mf](h)), ylim=c(min(c(sll,diff2)),max(c(sul,diff2))))
lines(r,sll,type='l',lty=2)
lines(r,sul,type='l',lty=2)


