# Michel Goulard's code
# notes de deperissement de lavande
fusaes <- as.matrix(read.table("lavande.txt"))

#on met les lignes de plantation en ligne
fusaes <- t(fusaes)

# on supprime la premiere ligne et les colonnes >=81
# pour avoir une plantation "rectangulaire" << Pedago>>
fusaes <- fusaes[-1,]
fusaes <- fusaes[,1:80]

# on regarde combien de chaque
#hist(fusaes,breaks=(0:6)-0.5,main="notes de dépérissement observées",xlab="frequence",ylab="notes")
hist(fusaes,breaks=(0:6)-0.5,main=" ",xlab="frequence",ylab="notes")


# on regarde ou c'est place
par(mar = c(3,3,1,1))

# we insert a blank row
new<-fusaes[1,]
new<-rbind(new,rep(NA,80)) 
for ( i in 2:11) {
	new<-rbind(new,fusaes[i,]) 
	new<-rbind(new,rep(NA,80)) 
}
library(fields)
 image.plot(1:80,1:22,t(new),xlab='x',ylab='y',col =gray((1:5)/5),axes=FALSE,legend.shrink=0.25)

# un variogramme le long des lignes
gam.fct <- function(tab,nd)
  {
    nx <- dim(tab)[1]
    ny <- dim(tab)[2]
    du <- rep(0,nd+1)
for (i in 1:nd)
  {
val1 <- tab[,(1+i):ny]
val2 <- tab[,1:(ny-i)]
du[i+1] <- sum((val1-val2)**2)
du[i+1] <- du[i+1]/(nx*(ny-i))
}
return(du)  }

#on l'essaie
par(mar = c(3,3,1,1))
plot(gam.fct(fusaes,20),type="l",main="essai variog log ligne")


# permutation complete
perm1.fct <- function(tab,gam.fct,n,nrep)
  {
    nx <- dim(tab)[1]
    ny <- dim(tab)[2]
vobs <- gam.fct(tab,n)
nz <- length(vobs)
tabres <- matrix(0,nz,nrep)
for (i in 1:nrep)
  {
tabres[,i] <- gam.fct(matrix(sample(tab),nx,ny),n)
  }
for (i in 1:nz)
  {
tabres[i,] <-sort(tabres[i,])
  }

 alphasur2 <-0.025
    k1 <- floor(alphasur2*nrep)
    k2 <- ceiling((1-alphasur2)*nrep)
zv <- range(c(vobs,tabres))
    plot(0:n,vobs,ylim=zv,xlab="distance",ylab="variogramme le long de la ligne",main=" ",pch=20)
    lines(0:n,vobs,lty=1)
    lines(0:n,tabres[,k1],lty=2)
    lines(0:n,tabres[,k2],lty=2)
    return()
  }


perm1.fct(fusaes,gam.fct,20,200)



# permutation des lignes
perm2.fct <- function(tab,gam.fct,n,nrep)
  {
    nx <- dim(tab)[1]
    ny <- dim(tab)[2]
vobs <- gam.fct(t(tab),n)
nz <- length(vobs)
tabres <- matrix(0,nz,nrep)
for (i in 1:nrep)
  {
tabres[,i] <- gam.fct(matrix(t(tab[sample(nx),]),nx,ny),n)
  }
for (i in 1:nz)
  {
tabres[i,] <-sort(tabres[i,])
  }

 alphasur2 <-0.025
    k1 <- floor(alphasur2*nrep)
    k2 <- ceiling((1-alphasur2)*nrep)
zv <- range(c(vobs,tabres))
    plot(0:n,vobs,ylim=zv,xlab="distance",ylab="variogramme entre lignes",main=" ",pch=20)
    lines(0:n,vobs,lwd=1)
    lines(0:n,tabres[,k1],lty=2)
    lines(0:n,tabres[,k2],lty=2)
    return()
  }
perm2.fct(fusaes,gam.fct,10,200)



# rotations des lignes
perm2b.fct <- function(tab,gam.fct,n,nrep)
  {
    nx <- dim(tab)[1]
    ny <- dim(tab)[2]
vobs <- gam.fct(t(tab),n)
nz <- length(vobs)
tabres <- matrix(0,nz,nrep)
for (i in 1:nrep)
  {
    t1 <- tab
    for (j in 1:nx)
      { id <- sample(ny,1)
        jd <- (1:ny)+id
        jd[jd>ny] <- jd[jd>ny]-ny
        t1[j,] <- t1[j,jd]
      }
tabres[,i] <- gam.fct(t(t1),n)
  }
for (i in 1:nz)
  {
tabres[i,] <-sort(tabres[i,])
  }

 alphasur2 <-0.025
    k1 <- floor(alphasur2*nrep)
    k2 <- ceiling((1-alphasur2)*nrep)
zv <- range(c(vobs,tabres))
    plot(0:n,vobs,ylim=zv,xlab="distance",ylab="variogramme entre lignes",main="rotations des lignes")
    lines(0:n,vobs,lwd=3)
    lines(0:n,tabres[,k1],lwd=2,lty=2)
    lines(0:n,tabres[,k2],lwd=2,lty=2)
    return()
  }

perm2b.fct(fusaes,gam.fct,10,200)


# permutation des seuls positifs (=attaques)
perm3.fct <- function(tab,gam.fct,n,nrep)
  {
    nx <- dim(tab)[1]
    ny <- dim(tab)[2]
vobs <- gam.fct(tab,n)
nz <- length(vobs)
tabres <- matrix(0,nz,nrep)
for (i in 1:nrep)
  {
    t1 <- tab
    t1[tab>0] <- sample(tab[tab>0])
tabres[,i] <- gam.fct(t1,n)
  }
for (i in 1:nz)
  {
tabres[i,] <-sort(tabres[i,])
  }

 alphasur2 <-0.025
    k1 <- floor(alphasur2*nrep)
    k2 <- ceiling((1-alphasur2)*nrep)
zv <- range(c(vobs,tabres))
    plot(vobs,ylim=zv,xlab="distance",ylab="variogramme le long de la ligne",main="randomization entre attaqués")
    lines(0:n,vobs,lwd=3)
    lines(0:n,tabres[,k1],lwd=2,lty=2)
    lines(0:n,tabres[,k2],lwd=2,lty=2)
    return()
  }

perm3.fct(fusaes,gam.fct,20,200)

# et si on permutait par rapport a la disposition dans les taches ?



par(mar = c(3,3,1,1))
perm4.fct <- function(tab,gam.fct,n,nrep)
  {
caldist.fct <- function(tab0)
{
  tab <- (-1)*tab0
  tab[tab==0] <- 1
  tab[tab<0] <- 0
  while(min(tab)==0){
  iu <- max(tab)
nx <- dim(tab)[1]
ny <- dim(tab)[2]
tt <- matrix(1,nx,ny+2)
tt[1:nx,2:(ny+1)] <- tab
ttm <- tt[,1:ny]
ttp <- tt[,3:(ny+2)]
ttpm <- ttm+ttp
sel <- (tab==0) & (ttpm>0)
  tab[sel] <- iu+1
}
  return(tab)
}

matdist <- caldist.fct(tab) -1
    nx <- dim(tab)[1]
    ny <- dim(tab)[2]
mdist <- max(matdist)
vobs <- gam.fct(tab,n)
nz <- length(vobs)
tabres <- matrix(0,nz,nrep)
for (i in 1:nrep)
  {
    t1 <- tab
    for (j in 1:mdist)
      {
        nu <- length(matdist[matdist==j])
        if(nu>1){
   t1[matdist==j] <- sample(tab[matdist==j])
   if(j <1){browser()}
   if(length(t1[matdist==j])-length(sample(tab[matdist==j]))){browser()}
      }}
    tabres[,i] <- gam.fct(t1,n)
  }
for (i in 1:nz)
  {
tabres[i,] <-sort(tabres[i,])
  }

 alphasur2 <-0.025
    k1 <- floor(alphasur2*nrep)
    k2 <- ceiling((1-alphasur2)*nrep)
zv <- range(c(vobs,tabres))
    plot(0:n,vobs,ylim=zv,xlab="distance",ylab="variogramme le long de la ligne",main="randomization entre attaqués, conditionnellement à la position")
    lines(0:n,vobs,lwd=3)
    lines(0:n,tabres[,k1],lwd=2,lty=2)
    lines(0:n,tabres[,k2],lwd=2,lty=2)
    return()
  }

perm4.fct(fusaes,gam.fct,20,200)


