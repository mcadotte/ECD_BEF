source("https://raw.githubusercontent.com/StefanoAllesina/Sao_Paulo_School/master/code/L-H.R")

Col <- c(rgb(13/255, 130/255, 97/255), rgb(173/255, 78/255, 23/255)) # Universal color combination

biomass <- function(A, K, p){
  B  <- matrix(0, nrow(p), ncol(p))
  for (i in 1:ncol(p)){
    pi       <- p[,i]
    Ai       <- -A[pi, pi, drop=FALSE]
    Ki       <- K[pi]
    B[pi, i] <- get_final_composition(Ai, Ki)
  }
  return(B)
}

n       <- 10 # Number of species
p       <- sapply(1:(2^n-1), function(k) as.logical(intToBits(k))[1:n]) # Combinations of species
p       <- p[,colSums(p)%in%c(1,2,4,8)]

quartz()
par(mar=c(2,2,1,1), mfrow=c(3,3))



# Control
K        <- rep(5, n)
A        <- 0.5 |> matrix(n, n)
diag(A)  <- 1
Control  <- biomass(A, K, p)
ContTot  <- sapply(c(1,2,4,8), function(s) Control[,colSums(p)==s] |> colSums() |> mean())
#ContCoef <- c(5, lm(I(ContTot-5) ~ log(c(1,2,4,8)) - 1)$coefficients)
ContCoef <- lm(ContTot ~ log(c(1,2,4,8)))$coefficients

# (A) ΔInter=0, ΔSlope<0
# ΔAlpha>0, Mean(ΔK)=0, Var(ΔK)=0
K          <- rep(5, n)
A          <- 0.69 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- lm(TreatTot ~ log(c(1,2,4,8)))$coefficients

plot(c(1,2,4,8),  ContTot, log="x", col=Col[1], pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") = 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") = 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") = 0")),pos=2)


plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") > 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") = 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") = 0")),pos=2)


# (B) ΔInter=0, ΔSlope>0
# ΔAlpha<0, Mean(ΔK)=0, Var(ΔK)=0
K          <- rep(5, n)
A          <- 0.35 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- lm(TreatTot ~ log(c(1,2,4,8)))$coefficients
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") < 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") = 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") = 0")),pos=2)

# (C) ΔInter<0, ΔSlope=0
# ΔAlpha=0, Mean(ΔK)<0, Var(ΔK)=medium
K          <- seq(1,5.5,l=10) # Mean(ΔK)<3.25 Var(ΔK)=2.29
A          <- 0.5 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- c(mean(K), lm(I(TreatTot-mean(K)) ~ log(c(1,2,4,8))-1)$coefficients)
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") = 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") < 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") > 0")),pos=2)

# (D) ΔInter<0, ΔSlope<0
# ΔAlpha=0, Mean(ΔK)<0, Var(ΔK)=small
K          <- rep(2, n) # Mean(ΔK)=1.5, Var(ΔK)=0
A          <- 0.5 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- c(mean(K), lm(I(TreatTot-mean(K)) ~ log(c(1,2,4,8))-1)$coefficients)
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(1,11,expression(paste(Delta,"mean(",alpha,") = 0")),pos=4)
text(1,10,expression(paste(Delta,"mean(",Kappa,") < 0")),pos=4)
text(1,9,expression(paste(Delta,"var(",Kappa,") = 0")),pos=4)

# (E) ΔInter<0, ΔSlope>0
# ΔAlpha=0, Mean(ΔK)<0, Var(ΔK)=large
K          <- exp(seq(-2.25, 2.25, l=10)) # Mean(ΔK)=2.40, Var(ΔK)=9.45
A          <- 0.5 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- lm(TreatTot ~ log(c(1,2,4,8)))$coefficients
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") = 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") < 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") >> 0")),pos=2)

# (F) ΔInter>0, ΔSlope=0
# ΔAlpha>0, Mean(ΔK)>0, Var(ΔK)=0
K          <- rep(7, n)
A          <- 0.59 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- lm(TreatTot ~ log(c(1,2,4,8)))$coefficients
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") > 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") > 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") = 0")),pos=2)

# (G) ΔInter>0, ΔSlope<0
# ΔAlpha>>0, Mean(ΔK)>0, Var(ΔK)=0
K          <- rep(7, n)
A          <- 0.75 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- lm(TreatTot ~ log(c(1,2,4,8)))$coefficients
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5,xlab="Number of species")
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") >> 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") > 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") = 0")),pos=2)

# (H) ΔInter>0, ΔSlope>0
# ΔAlpha=0, Mean(ΔK)>0, Var(ΔK)=0
K          <- rep(7, n)
A          <- 0.5 |> matrix(n, n)
diag(A)    <- 1
Treat      <- biomass(A, K, p)
TreatTot   <- sapply(c(1,2,4,8), function(s) Treat[,colSums(p)==s] |> colSums() |> mean())
TreatCoef  <- lm(TreatTot ~ log(c(1,2,4,8)))$coefficients
plot(rep(c(1,2,4,8),2),  c(ContTot, TreatTot), log="x", col=rep(Col,e=4), pch=19, ylim=c(2,13), bty="n", las=1, cex=1.5)
lines(c(1,8), c(ContCoef[1],ContCoef[1]+ContCoef[2]*log(8)), col=Col[1], lty=2, lwd=2.5)
lines(c(1,8), c(TreatCoef[1],TreatCoef[1]+TreatCoef[2]*log(8)), col=Col[2], lwd=2.5)
text(8,5,expression(paste(Delta,"mean(",alpha,") = 0")),pos=2)
text(8,4,expression(paste(Delta,"mean(",Kappa,") > 0")),pos=2)
text(8,3,expression(paste(Delta,"var(",Kappa,") = 0")),pos=2)



