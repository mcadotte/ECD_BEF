###using Chesson modern coexistence theory to assess ECD impacts on coexistence and inferences about BEF relationships. 

rho<-seq(0.25,0.99,0.01)
inv.rho<-1/rho
FR<-seq(0,4.2,0.01)

min<-NULL
max<-NULL
for (i in 1:length(rho)){
  tmp<-FR[FR>rho[i] & FR<inv.rho[i]]
  min[i]<-min(tmp)
  max[i]<-max(tmp)
}

ND<-1-rho

ND_FD_1<-function(r1,r2,a12,a21,a11,a22){
  rND<-sqrt((a12*a21)/(a11*a22))
  rFD<-(log(r1)/log(r2))*sqrt((a12*a11)/(a21*a22))
  return(data.frame(Niche_diff=rND,Fitness_ratio=rFD))
}

##this is a formulation including Ks if estimation of fitness and niche differences
ND_FD_2<-function(K1,K2,a12,a21){
  t.ND<-sqrt((a12/K1)*(a21/K2))
  t.FD<-(K1/K2)*sqrt((a21*K1)/(a12*K2))
  return(data.frame(Niche_diff=t.ND,Fitness_ratio=t.FD))
}

###using formulation with K

quartz()
plot(ND,max,type="l",ylim=c(0,4),lwd=3,
     xlab=expression("Niche difference (1-"*rho*")"),
     ylab=expression("Fitness ratio ("*kappa*"1/"*kappa*"2)"))
polygon(c(min(ND),max(ND),max(ND),ND),c(1,min(min),max(max),max),col="grey85")

lines(ND,max,type="l",lwd=3)
lines(ND,min,type="l",lwd=3)

K1<-15; K2<-10; a12<-3; a21<-2

p1<-ND_FD_2(K1,K2,a12,a21)
points(p1,pch=19,cex=1.6,col="darkslategray4")

K1<-12
p2<-ND_FD_2(K1,K2,a12,a21)
points(p2,pch=19,cex=1.6,col="darkorange3")

arrows(p1[,1],p1[,2],p2[,1],p2[,2],lwd=2,length=0.15,angle=20)

##next case
K1<-12; K2<-8; a12<-6; a21<-2
p3<-ND_FD_2(K1,K2,a12,a21)
p3

points(p3,pch=19,cex=1.6,col="darkslategray4")

K2<-5.3
p4<-ND_FD_2(K1,K2,a12,a21)
p4
points(p4,pch=19,cex=1.6,col="darkorange3")
arrows(p3[,1],p3[,2],p4[,1],p4[,2],lwd=2,length=0.15,angle=20)


text(p1,"A",pos=4)
text(p4[1],p4[2]+0.2,"B")
text(0.5,1, pos=4, "Coexistence",cex=1)
text(0.03,0, pos=4, "Competative exclusion",cex=1)
text(0.2,3, pos=4, "Competative exclusion",cex=1)


###ND
K1<-14; K2<-10; a12<-6; a21<-9

N1<-ND_FD_2(K1,K2,a12,a21)
points(N1,pch=19,cex=1.6,col="darkslategray4")

a12<-4.7;a21<-6.7
N2<-ND_FD_2(K1,K2,a12,a21)
N2
points(N2,pch=19,cex=1.6,col="darkorange3")
arrows(N1[,1],N1[,2],N2[,1],N2[,2],lwd=2,length=0.15,angle=20)
text(N1[,1],(N1[,2]+0.1),"C",pos=4)


K1<-14; K2<-8; a12<-5; a21<-6

N3<-ND_FD_2(K1,K2,a12,a21)
N3
points(N3,pch=19,cex=1.6,col="darkslategray4")

a12<-6; a21<-9
N4<-ND_FD_2(K1,K2,a12,a21)
N4

points(N4,pch=19,cex=1.6,col="darkorange3")
arrows(N3[,1],N3[,2],N4[,1],N4[,2],lwd=2,length=0.15,angle=20)
text(N4[,1],(N4[,2]+0.1),"D",pos=4)



###standard version without K in estimating ND and FD
quartz()
plot(ND,max,type="l",ylim=c(0,4),lwd=3,
     xlab=expression("Niche difference (1-"*rho*")"),
     ylab=expression("Fitness ratio ("*kappa*"1/"*kappa*"2)"))
polygon(c(min(ND),max(ND),max(ND),ND),c(1,min(min),max(max),max),col="grey85")

lines(ND,max,type="l",lwd=3)
lines(ND,min,type="l",lwd=3)

r1<-1.4; r2<-1.2; a11<-0.6; a22<-0.6; a12<-0.3; a21<-0.2

p1<-ND_FD(r1,r2,a12,a21,a11,a22)
points(p1,pch=19,cex=1.6,col="darkslategray4")

r1<-1.2
p2<-ND_FD(r1,r2,a12,a21,a11,a22)
points(p2,pch=19,cex=1.6,col="darkorange3")

arrows(p1[,1],p1[,2],p2[,1],p2[,2],lwd=2,length=0.15,angle=20)

r1<-1.3; a12<-0.2;a21<-0.45;a11<-0.9;a22<-0.9
p3<-ND_FD(r1,r2,a12,a21,a11,a22)
p3

points(p3,pch=19,cex=1.6,col="darkslategray4")

r2<-1.6
p4<-ND_FD(r1,r2,a12,a21,a11,a22)
p4
points(p4,pch=19,cex=1.6,col="darkorange3")
arrows(p3[,1],p3[,2],p4[,1],p4[,2],lwd=2,length=0.15,angle=20)


text(p1,"1.3",pos=4)
text(p4,"1.4",pos=4)
text(0.5,1, pos=4, "Coexistence",cex=0.8)
text(0.03,0, pos=4, "Competative exclusion",cex=0.8)
text(0.2,3, pos=4, "Competative exclusion",cex=0.8)


###ND
r1<-1.45; r2<-1.2; a11<-0.45; a22<-0.45; a12<-0.3; a21<-0.2

N1<-ND_FD(r1,r2,a12,a21,a11,a22)
N1

a21<-0.25
N2<-ND_FD(r1,r2,a12,a21,a11,a22)
N2
points(N1,pch=19,cex=1.6,col="darkslategray4")
points(N2,pch=19,cex=1.6,col="darkorange3")
arrows(N1[,1],N1[,2],N2[,1],N2[,2],lwd=2,length=0.15,angle=20)
text(N1[,1],(N1[,2]+0.1),"3.1",pos=4)


r1<-1.55; r2<-1.2; a11<-0.35; a22<-0.35; a12<-0.3; a21<-0.2

N3<-ND_FD(r1,r2,a12,a21,a11,a22)
N3

a21<-0.15
N4<-ND_FD(r1,r2,a12,a21,a11,a22)
N4
points(N3,pch=19,cex=1.6,col="darkslategray4")
points(N4,pch=19,cex=1.6,col="darkorange3")
arrows(N3[,1],N3[,2],N4[,1],N4[,2],lwd=2,length=0.15,angle=20)
text(N4[,1],(N4[,2]+0.1),"3.1",pos=4)
