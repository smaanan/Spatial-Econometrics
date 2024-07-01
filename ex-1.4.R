library(RandomFields)
par(pty='s')
cols<- gray(0:16 / 16)
set.seed(19)
 model <- "stable"   
      mean <- 0
      variance <- 4
      nugget <- 1
      scale <- 10
      alpha <- 1  
      step <- 1   
      x <- seq(0, 10, length=100) 
      y <- seq(0, 10, length=100)     
  RFparameters(TBM2.lines=2)
      f <- GaussRF(method="TBM3",x=x, y=y, model=model, grid=TRUE,param=c(mean, variance, nugget, scale, alpha))
		         image(x, y, f,col=cols)


set.seed(19)
 model <- "stable"   
      mean <- 0
      variance <- 4
      nugget <- 0
      scale <- 10
      alpha <- 1  
      step <- 1   
      x <- seq(0, 10, length=100) 
      y <- seq(0, 10, length=100)     
  RFparameters(TBM2.lines=10)
set.seed(19)
      f <- GaussRF(method="TBM3",x=x, y=y, model=model, grid=TRUE, 
                   param=c(mean, variance, nugget, scale, alpha))
par(pty='s')
		         image(x, y, f,col=cols)


library(RandomFields)	
set.seed(1)
param <- c(0, 1, 0, 1)
model <- "exponential"
#RFparameters(PracticalRange=FALSE)
krige.method <- "O" ## random field assumption corresponding to
                    ## those of ordinary kriging
length<-100
x <-  seq(0, 1, length=length)
y <-  seq(0, 1, length=length)
xy<-expand.grid(x,y)
data <- GaussRF(x=xy[,1], y=xy[,2], grid=FALSE, model=model, param=param)
# another grid, where values are to be simulated
# standardisation of the output
lim <- range( x)
zlim <- c(min(data)-0.1, max(data)+0.1)
colour <- gray( seq(0,1,l= 16))
data.mat<-matrix(data,length(x),length(x))
## visualise generated spatial data
par(pty='s')
image(x, y, data.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)



set.seed(19)
par(pty='s')
sub.pos<-sample(1:length(data),100)
plot(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2],xlab='x',ylab='y',type='n')
points(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2],pch=20)



kdata <-  Kriging(krige.method=krige.method,
              x=xy[,1], y=xy[,2],  grid=FALSE,
              model=model, param=param,
              given=xy[sub.pos[1:25],], data=data[sub.pos[1:25]])
kdata.mat<-matrix(kdata,length(x),length(x))
par(pty='s')

image(x, y, kdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2])




#conditional simulation n=25
cdata <- CondSimu(krige.method, x=xy[,1], y=xy[,2],  grid=FALSE,
               model=model, param=param,
               given=xy[sub.pos[1:25],],data=data[sub.pos[1:25]])

cdata.mat<-matrix(cdata,length(x),length(x))
par(pty='s')

image(x, y, cdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:25],1],xy[sub.pos[1:25],2])



#conditional simulation n=50

cdata <- CondSimu(krige.method, x=xy[,1], y=xy[,2],  grid=FALSE,
               model=model, param=param,given=xy[sub.pos[1:50],],data=data[sub.pos[1:50]])


cdata.mat<-matrix(cdata,length(x),length(x))
par(pty='s')

image(x, y, cdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:50],1],xy[sub.pos[1:50],2])



#conditional simulation n=100
cdata <- CondSimu(krige.method, x=xy[,1], y=xy[,2],  grid=FALSE,
               model=model, param=param,given=xy[sub.pos[1:100],],data=data[sub.pos[1:100]])

cdata.mat<-matrix(cdata,length(x),length(x))
par(pty='s')

image(x, y, cdata.mat, xlim=lim, ylim=lim, zlim=zlim, col=colour)
points(xy[sub.pos[1:100],1],xy[sub.pos[1:100],2])


