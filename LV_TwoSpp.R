###two species LV model

#function for two-species LV discrete time steps
LV2<-function(r1,r2,N1,N2,a12,a21,K1,K2,G=100){
  
  out.N1<-NULL
  out.N2<-NULL
  out.N1[1]<-N1
  out.N2[1]<-N2
  Total<-N1+N2
  
  for (i in 1:G){
    out.N1[i+1]<-
      out.N1[i]+(r1*out.N1[i]*(1-((out.N1[i]+a12*out.N2[i])/K1)))
    out.N2[i+1]<-
      out.N2[i]+(r2*out.N2[i]*(1-((out.N2[i]+a21*out.N1[i])/K2)))
    Total[i+1]<-out.N1[i+1]+out.N2[i+1]
  }
  return(as.data.frame(cbind(out.N1,out.N2, Total)))
}


###ECD effects on r -lower
r1<-0.2
r2<-0.18
N1<-5
N2<-5
a12<-0.1
a21<-0.1
K1<-50
K2<-60


r.out1<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

quartz()
par(mfrow=c(1,3))
plot(1:nrow(r.out1),r.out1$Total,ylim=c(0,max(r.out1$Total)),
     type="l", lwd=3,col="grey70", xlab="Time",ylab="Abundance",
     main = "a) ECD decreases r",cex.lab=1.5,cex.main=1.5)
points(1:nrow(r.out1),r.out1$out.N1,type="l", lwd=3,col="grey70")
points(1:nrow(r.out1),r.out1$out.N2,type="l", lwd=3,col="grey70")

r1<-0.1
r2<-0.09
r.out1.L<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

points(1:nrow(r.out1.L),r.out1.L$Total,type="l", lwd=3,col="darkslategray")
points(1:nrow(r.out1.L),r.out1.L$out.N1,type="l", lwd=3,col="darkslategray4")
points(1:nrow(r.out1.L),r.out1.L$out.N2,type="l", lwd=3,col="darkslategray3")

text(70,95,pos=4,"Total",col="darkslategray")
text(70,60,pos=4,"Species 2",col="darkslategray3")
text(70,40,pos=4,"Species 1",col="darkslategray4")


x1<-min((1:length(r.out1$Total))[r.out1$Total>=0.95*max(r.out1$Total)])
x2<-min((1:length(r.out1.L$Total))[r.out1.L$Total>=0.95*max(r.out1.L$Total)])

polygon(c(x1,x1,x2,x2),c(0,2,2,0),col="black")


###ECD effects on r -increase
r1<-3
r2<-2

r.out1.H<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

plot(1:nrow(r.out1),r.out1$Total,
     type="l", lwd=3,col="grey70", xlab="Time",ylab="Abundance",
     main = "b) ECD increases r",ylim=c(0,130),cex.lab=1.5,cex.main=1.5)
points(1:nrow(r.out1),r.out1$out.N1,type="l", lwd=3,col="grey70")
points(1:nrow(r.out1),r.out1$out.N2,type="l", lwd=3,col="grey70")


points(1:nrow(r.out1.H),r.out1.H$Total,type="l", lwd=3,col="darkslategray")
points(1:nrow(r.out1.H),r.out1.H$out.N1,type="l", lwd=3,col="darkslategray4")
points(1:nrow(r.out1.H),r.out1.H$out.N2,type="l", lwd=3,col="darkslategray3")

#do plot r by CV =sd/mean from gen 20 on. Kept ratio of r constant, r2 = r1*0.8

CV<-function(x,na.rm=TRUE){
  return(sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))
}

rs<-seq(0.3,3,0.01)
CVs<-NULL
for(i in 1:length(rs)){
  tmp<-LV2(rs[i],rs[i]*0.8,N1,N2,a12,a21,K1,K2)
  tmp<-na.omit(tmp)
  tmp<-tmp[41:nrow(tmp),]
  CVs[i]<-CV(tmp$Total)
}
CVs2<-CVs[CVs!=-Inf]

plot(rs[CVs!=-Inf],CVs2, type="l", lwd=3,col="darkslategray",
     ylab="Coefficient of variation",
     xlab="Intrinsic rate of increase",
     main="c) Coefficient of variation",cex.lab=1.5,cex.main=1.5)



  
#examine ECD effects on competition 

r1<-0.2
r2<-0.18
N1<-5
N2<-5
a12<-0.2
a21<-0.2
K1<-50
K2<-60

bl<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

#increases in competition
a12<-0.05
a12<-0.03
CO1<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

#quartz()
#par(mfrow=c(1,2))
plot(1:nrow(bl),bl$Total,ylim=c(0,max(KO1$Total)),
     type="l", lwd=3,col="grey70", xlab="Time",ylab="Abundance",
     main = "GCD changes the stength of competition")
points(1:nrow(bl),bl$out.N1,type="l", lwd=3,col="grey70")
points(1:nrow(bl),bl$out.N2,type="l", lwd=3,col="grey70")

points(1:nrow(CO1),CO1$Total,type="l", lwd=3,col="darkslategray")
points(1:nrow(CO1),CO1$out.N1,type="l", lwd=3,col="darkslategray4")
points(1:nrow(CO1),CO1$out.N2,type="l", lwd=3,col="darkslategray3")


#Average C decreases
a12<-0.65
a21<-0.70
CO2<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

points(1:nrow(CO2),CO2$Total,type="l", lwd=3,col="darkorange4")
points(1:nrow(CO2),CO2$out.N1,type="l", lwd=3,col="darkorange3")
points(1:nrow(CO2),CO2$out.N2,type="l", lwd=3,col="darkorange2")


#increasing asymmetry 
a12<-0.05
a21<-0.35
CO3<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

plot(1:nrow(bl),bl$Total,ylim=c(0,max(KO1$Total)),
     type="l", lwd=3,col="grey70", xlab="Time",ylab="Abundance",
     main = "a) GCD changes the symmetry of competition")
points(1:nrow(bl),bl$out.N1,type="l", lwd=3,col="grey70")
points(1:nrow(bl),bl$out.N2,type="l", lwd=3,col="grey70")

points(1:nrow(CO3),CO3$Total,type="l", lwd=3,col="darkslategray")
points(1:nrow(CO3),CO3$out.N1,type="l", lwd=3,col="darkslategray4")
points(1:nrow(CO3),CO3$out.N2,type="l", lwd=3,col="darkslategray3")




###ECD effects on carrying capacity

r1<-0.2
r2<-0.18
N1<-5
N2<-5
a12<-0.1
a21<-0.1
K1<-50
K2<-60

bl<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

#Average K increases
K1<-60
K2<-70
KO1<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

#quartz()
#par(mfrow=c(1,2))
#plot(1:nrow(bl),bl$Total,ylim=c(0,max(KO1$Total)),
#     type="l", lwd=3,col="grey70", xlab="Time",ylab="Abundance",
#     main = "a) GCD changes average K")
#points(1:nrow(bl),bl$out.N1,type="l", lwd=3,col="grey70")
#points(1:nrow(bl),bl$out.N2,type="l", lwd=3,col="grey70")

#points(1:nrow(KO1),KO1$Total,type="l", lwd=3,col="darkslategray")
#points(1:nrow(KO1),KO1$out.N1,type="l", lwd=3,col="darkslategray4")
#points(1:nrow(KO1),KO1$out.N2,type="l", lwd=3,col="darkslategray3")


#Average K decreases
K1<-40
K2<-50
KO2<-LV2(r1,r2,N1,N2,a12,a21,K1,K2)

#points(1:nrow(KO2),KO2$Total,type="l", lwd=3,col="darkorange4")
#points(1:nrow(KO2),KO2$out.N1,type="l", lwd=3,col="darkorange3")
#points(1:nrow(KO2),KO2$out.N2,type="l", lwd=3,col="darkorange2")

#K variance increases no average effect
#not shown

#K variance increases decrease in average
#simulate across range of K values will average = 45 
#and high variance then heat map showing magnitude of effect
#on community biomass as tmp$total - bl$total


r1<-0.2
r2<-0.18
N1<-5
N2<-5
a12<-0.1
a21<-0.1

tot<-data.frame(Mean=0,SD=0,Slope=0,Int=0)
SD<-1:15
Mu<-45

for (j in 1:length(SD)){
  slope.tmp<-NULL
  int.tmp<-NULL
  for (i in 1:99){
    tmp<-rnorm(2,Mu,SD[j])
    t.mod<-LV2(r1,r2,N1,N2,a12,a21,tmp[1],tmp[2])
    t.an<-data.frame(SR=c(1,1,2),BM=abs(c(t.mod$out.N1[101],t.mod$out.N2[101],t.mod$Total[101])))
    t.lm<-lm(BM~SR,data=t.an)
    slope.tmp[i]<-summary(t.lm)$coefficients[2,1]
    int.tmp[i]<-summary(t.lm)$coefficients[1,1]
  }
  
  tot[j,]<-c(Mu,SD[j],mean(slope.tmp),mean(int.tmp))
} 

  