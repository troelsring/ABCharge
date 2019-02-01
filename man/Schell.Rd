\name{Schell}
\alias{Schell}
\alias{S1}
\alias{Phosphate}
\alias{Citrate}
\docType{data}
\title{
Data from Schell et al 2018
}
\description{
Data from Schell et al J Cem End Data 2018; 63: 2151-2161
}
\usage{data("Schell")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{Na,K,Cl,Ca,Mg,Lact}}{a numeric vector of strong ions}
    \item{\code{Alb,PCO2}}{a numeric vector of Albumin and PCO2, 0 if absent}
    \item{\code{WA}}{list of buffers and pKas}
    \item{\code{TOTAL}}{Matrix of buffer concentrations}
    \item{\code{S1}}{Results from Schell et al}
    \item{\code{A,b}}{Constants from Davies' equation}
  }
}
\details{
A series of 30 mixtures of citric acid with dihydrogen sodium phosphate
}
\source{Schell et al. J Cem End Data 2018; 63: 2151-2161

}

\examples{

data(Schell)
pHRES <- c();KpH <- c(); Charge <- IonStr <- iterations <- CCCH <-
UNCOR_pH <-  c()
FF1 <- array(rep(NA,3*30),c(30,3))
SPECS <- array(rep(NA,8*30),c(30,8))

for (j in 1:30) {

  UNCOR_pH[j] <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$UNCOR_pH
  pHRES[j] <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$pH
  SPECS[j,] <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$SPECS
  FF1[j,]  <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$FF1
  Charge[j]<- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$Charge
 }
 ###**Point 1**

#### Analysis in ideal state

pHOBS <- S1[,5]


DFF <- data.frame(pHOBS,UNCOR_pH)
PP <- with(DFF,cor.test(pHOBS,UNCOR_pH)) #0.9999631
PP

gg <- ggplot2::ggplot()+ggplot2::geom_point(data=DFF,ggplot2::aes(
x=pHOBS,y=UNCOR_pH),
size=3) + ggplot2::geom_abline(intercept=0,slope=1,col="red")+
  ggplot2::scale_x_continuous(limits=c(2,9))+ggplot2::scale_y_continuous(
  limits=c(2,9))+ggplot2::xlab("pH reported by Schell")+
  ggplot2::ylab("pH calculated from charge-balance")+
  ggplot2::annotate("text",x=4,y=7,label="Correlation = 0.9999631",size=7.5)
gg + ggplot2::ggtitle("pH from Schell, S1 column 5 on abscissa.Activity
coefficients forced = 1")+
  ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold"))

#options(digits=10)
ss <-pHOBS-UNCOR_pH ;
ss <- summary(ss)

###**Point 2**

#### Analysis with correction for ionic strength

pHOBS = S1[,7]
cor.test(pHOBS,pHRES) #0.999959

DFF <- data.frame(pHOBS,pHRES)
gg <- ggplot2::ggplot()+ggplot2::geom_point(data=DFF,ggplot2::aes(
x=pHOBS,y=pHRES),size=3) + ggplot2::geom_abline(intercept=0,slope=1,col="red")+
  ggplot2::scale_x_continuous(limits=c(2,8))+ggplot2::scale_y_continuous(
  limits=c(2,8))+ggplot2::xlab("pH reported by Schell")+
  ggplot2::ylab("pH calculated from charge-balance")
gg + ggplot2::annotate("text",5,7,label="Correlation = 0.99996",size=7.5)+
  ggplot2::ggtitle("pH from charge-balance and
          as reported (S1, column 7) by Schell")+
  ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold"))

###**Point 3**
####Analysis of pKs from Schell et al^1^

IS_stated <- S1[,13]

A <- 0.5108
f1  <- 10^(-A*1^2*(sqrt(IS_stated)/(1+sqrt(IS_stated))-IS_stated*0.1))
f2  <- 10^(-A*2^2*(sqrt(IS_stated)/(1+sqrt(IS_stated))-IS_stated*0.1))
f3  <- 10^(-A*3^2*(sqrt(IS_stated)/(1+sqrt(IS_stated))-IS_stated*0.1))

IND <- Citrate[,6]>0 & Citrate[,3]>0;   A1 <- sum(IND)
CitpK1 <- (S1[,8]-log10(Citrate[,6]*f1/Citrate[,3]))[IND]
IS1 <- IS_stated[IND]
summary(CitpK1)

IND <- Citrate[,6]>0 & Citrate[,9]>0;A2 <- sum(IND)
CitpK2 <- (S1[,8]-log10(Citrate[,9]*f2/(Citrate[,6]))*f1)[IND]
IS2 <- IS_stated[IND]
summary(CitpK2)

IND <- Citrate[,12]>0 & Citrate[,9]>0;A3 <- sum(IND)
CitpK3 <- (S1[,8]-log10(Citrate[,12]*f3/(Citrate[,9]))*f2)[IND]
IS3 <- IS_stated[IND]
summary(CitpK3)

CIT <- c(CitpK1,CitpK2,CitpK3)
ISS <- c(IS1,IS2,IS3)
PK <- factor(c(rep(1,A1),rep(2,A2),rep(3,A3)),labels=c("pK1","pK2","pK3"))
PKT <- c(rep(3.13,A1),rep(4.76,A2),rep(6.4,A3))

DDD <- data.frame(CIT,PK,ISS,PKT)

ff <- ggplot2::ggplot(DDD,ggplot2::aes(PK,CIT))+ggplot2::geom_boxplot() +
ggplot2::geom_dotplot(binaxis="y",stackdir="center",binwidth=.1)+
ggplot2::xlab("")+
ggplot2::ylab("Citrate pKa1-3")
ff + ggplot2::geom_segment(x=0.75,xend=1.25,y=3.13,yend=3.13,col="red",size=2)+
  ggplot2::geom_segment(x=1.75,xend=2.25,y=4.76,yend=4.76,col="red",size=2)+
  ggplot2::geom_segment(x=2.75,xend=3.25,y=6.40,yend=6.40,col="red",size=2)+
  ggplot2::ggtitle("The 3 observed Citrate pK values based on\n
  measured H activity \n and printed species distribution.
  3.13, 4.76, and 6.40 indicated by red bars")
#ggplot2::ggsave("Citrate-pKs.pdf",width = 8.3, height = 8.3, units = c("cm"))

dd <- ggplot2::ggplot()
dd <- dd + ggplot2::geom_point(data=DDD,ggplot2::aes(x=ISS,y=PKT-CIT,col=PK),
size=2)+ggplot2::xlab("Ionic strength")+ggplot2::ylab("True- calculated pK")
dd



####Phosphate


IND <- Phosphate[,6]>0 & Phosphate[,3]>0; B1 <- sum(IND)
PhospK1 <- (S1[,8]-log10(Phosphate[,6]*f1/Phosphate[,3]))[IND]
summary(PhospK1)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 2.020   2.142   2.158   2.162   2.169   2.409
#
IS1 <- IS_stated[IND]

IND <- Phosphate[,6]>0 & Phosphate[,9]>0; B2 <- sum(IND)
PhospK2 <- (S1[,8]-log10(Phosphate[,9]*f2/(Phosphate[,6]))*f1)[IND]
summary(PhospK2)
  # Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 6.313   6.610   6.913   6.942   7.232   7.613
IS2 <- IS_stated[IND]
IND <- Phosphate[,12]>0 & Phosphate[,9]>0; B3 <- sum(IND)
PhospK3 <- (S1[,8]-log10(Phosphate[,12]*f3/(Phosphate[,9]))*f2)[IND]
summary(PhospK3)
  #  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
  # 8.848   8.848   8.848   8.848   8.848   8.848
IS3 <- IS_stated[IND]

  PHO <- c(PhospK1,PhospK2,PhospK3)
PK <- factor(c(rep(1,B1),rep(2,B2),rep(3,B3)),labels=c("pK1","pK2","pK3"))
PKT <- c(rep(2.16,B1),rep(7.21,B2),rep(12.32,B3))
ISS <- c(IS1,IS2,IS3)
DDD <- data.frame(PHO,PK,PKT,ISS)

ff <- ggplot2::ggplot(DDD,ggplot2::aes(PK,PHO))+ggplot2::geom_boxplot() +
ggplot2::ylim(c(1,13))

ff + ggplot2::geom_segment(x=0.75,xend=1.25,y=2.16,yend=2.16,
col="blue",size=2)+
  ggplot2::geom_segment(x=1.75,xend=2.25,y=7.21,yend=7.21,col="blue",size=2)+
  ggplot2::geom_segment(x=2.75,xend=3.25,y=12.32,yend=12.32,col="blue",size=2)+
  ggplot2::geom_dotplot(binaxis="y",stackdir="center",binwidth=.2)+
  ggplot2::xlab("")+ggplot2::ylab("Phosphate pKa1-3")+
  ggplot2::ggtitle("The 3 observed Phosphate pK values based on measured
  H activity and printed species distribution. 2.16, 7.21, and 12.32
  indicated by blue bars")

####Charge balance
Charge <- S1[,3]*2-Citrate[,6]*1e-3-Citrate[,9]*2*1e-3-Citrate[,12]*3*
1e-3-Phosphate[,6]*1e-3-Phosphate[,9]*2*1e-3-Phosphate[,12]*3*1e-3 + 10^-S1[,7]-
kw/10^-S1[,7]
(ss <- summary(Charge))
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0009095 -0.0002690  0.0001403  0.0002145  0.0008153  0.0011996
GRP <- factor(rep(1,length(Charge)))
DFF <- data.frame(Charge,GRP)
ff <- ggplot2::ggplot(DFF,ggplot2::aes(GRP,Charge))+ggplot2::geom_boxplot() +
ggplot2::geom_dotplot(binaxis="y",stackdir="center",binwidth=.00005)+
ggplot2::xlab("")+ggplot2::ylab("Charge Eql/l")
ff

##Our charge balance
ZZ <- c(0,-1,-2,-3,0,-1,-2,-3)

CH <- c()

for (i in 1:30) CH[i] <- sum(SPECS[i,] * ZZ) + Na[i] + 10^-pHRES[i] - kw/(10^-pHRES[i]*FF1[i,1]^2)
(ss <- summary(CH))

## Revocery of pks
# Is total phosphate conserved
PHOSTOT <- c()
for (i in 1:30) PHOSTOT[i] <- sum(SPECS[i,1:4])

 all.equal(PHOSTOT,S1[,3]) #TRUE
# Is total citrate conserved
CITTOT <- c()
for (i in 1:30) CITTOT[i] <- sum(SPECS[i,5:8])

all.equal(CITTOT,S1[,2]) #TRUE


RPK1 <- (-log10(SPECS[,2]*FF1[,1]^2*10^-pHRES/SPECS[,1]))
RPK2 <- (-log10(SPECS[,3]*FF1[,2]*10^-pHRES/SPECS[,2]))
RPK3 <- (-log10(SPECS[,4]*FF1[,1]*FF1[,3]*10^-pHRES/(SPECS[,3]*FF1[,2])))

RCK1 <- (-log10(SPECS[,6]*FF1[,1]^2*10^-pHRES/SPECS[,5]))
RCK2 <- (-log10(SPECS[,7]*FF1[,2]*10^-pHRES/SPECS[,6]))
RCK3 <- (-log10(SPECS[,8]*FF1[,1]*FF1[,3]*10^-pHRES/(SPECS[,7]*FF1[,2])))

summary(RPK1)
summary(RPK2)
summary(RPK3)
summary(RCK1)
summary(RCK2)
summary(RCK3)

}

\keyword{datasets}
