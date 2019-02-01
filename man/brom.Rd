\name{brom}
\alias{brom}
\alias{TOTAL}
\docType{data}
\title{
bromic acid data from Glaser et al 2014
}
\description{
  Data to find optimal pK for bromic acid
}
\usage{data("brom")}
\format{
 A workspace.
  \describe{
    \item{brom}{Dataframe of Table S4 in Glaser et al.}
    \item{TOTAL}{Matrix of buffer concentrations.}
  }
}
\source{
Glaser RE et al. Dynamical approach to multiequilibria problems for mixtures of
 acids and their conjugate bases. J Chem Educ 2014; 91: 1009-1016
}
\examples{
data(brom)
pK11 <- -3
pK12 <- 2
Na <- 0.1
kw <- 1e-14
PKBr <- seq(-2,2,by=0.2)
KABr <- 10^-PKBr

pHRES <- c();KpH <- c(); Charge <- IonStr <- iterations <- CCCH <- RCCH <-
    UNCOR_pH <- III <- c()
FF1 <- array(rep(NA,2*21),c(21,2))
SPECS <- FSPECS <- array(rep(NA,5*21),c(21,5))

for (j in 1:21) {
    WA <- list(buffs=list(),KA=list())
    WA$buffs = list("Sulfuric","Brom")
    WA$KA <- list(c(10^-pK11,10^-pK12),c(KABr[j]))

    pHRES[j]  <-  pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
                             TOT=TOTAL[j,],A=0.5108,b=0.1)$pH
    (UNCOR_pH[j] <- pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
    WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$UNCOR_pH)
    SPECS[j,] <- pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
                            TOT=TOTAL[j,],A=0.5108,b=0.1)$SPECS
    FF1[j,]  <- pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
                           TOT=TOTAL[j,],A=0.5108,b=0.1)$FF1
    Charge[j]<- pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
                           TOT=TOTAL[j,],A=0.5108,b=0.1)$Charge
    III[j]  <- pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
                          TOT=TOTAL[j,],A=0.5108,b=0.1)$II
    FSPECS[j,] <- pH_general(Na=Na,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
                             TOT=TOTAL[j,],A=0.5108,b=0.1)$FSPECS
}

ddf <- data.frame(pHobs=brom$pH,pHmod=UNCOR_pH)
gg <- ggplot2::ggplot()
gg <- gg + ggplot2::geom_point(data=ddf,ggplot2::aes(x=pHobs,y=pHmod))
gg + ggplot2::geom_abline(intercept=0,slope=1,col="red")+
ggplot2::xlab("pH from Glaser")+
    ggplot2::ylab("Modeled pH")

#correlation modeled-measured
(XX <- cor.test(brom$pH,UNCOR_pH)$estimate)

summary(brom$pH-UNCOR_pH)
GLBromat <- brom[,4]
GLBrH    <- brom[,5]

dddd <- data.frame(pH=UNCOR_pH,PKBR=PKBr,BRO3=FSPECS[,5],HBRO3=FSPECS[,4],
                   GLBromat,GLBrH)

gg <- ggplot2::ggplot()
gg <- gg + ggplot2::geom_line(data=dddd,ggplot2::aes(x=PKBr,y=pH),col="green",
                              size=1)
gg <- gg + ggplot2::geom_line(data=dddd,ggplot2::aes(x=PKBr,y=BRO3),col="red",
                              size=1)
gg <- gg + ggplot2::geom_line(data=dddd,ggplot2::aes(x=PKBr,y=HBRO3),col="blue",
                              size=1)
gg <- gg + ggplot2::geom_point(data=dddd,ggplot2::aes(x=PKBr,y=pH),col="green",
                               size=2)
gg <- gg + ggplot2::geom_point(data=dddd,ggplot2::aes(x=PKBr,y=BRO3),col="red",
                               size=2)
gg <- gg + ggplot2::geom_point(data=dddd,ggplot2::aes(x=PKBr,y=HBRO3),
                               col="blue",size=2)+
    ggplot2::scale_y_continuous(breaks=seq(-0.02,0.12,by=0.02))+
    ggplot2::scale_x_continuous(breaks=seq(-2.5,2.5,by=0.5))+
    ggplot2::xlab("pKBr")+
    ggplot2::ylab("pH (green) and concentrations M")
gg

dd <- ggplot2::ggplot()
dd <- dd + ggplot2::geom_point(data=dddd,ggplot2::aes(x=GLBromat,y=BRO3))+
    ggplot2::geom_abline(intercept=0,slope=1,col="red")
dd

dd <- ggplot2::ggplot()
dd <- dd + ggplot2::geom_point(data=dddd,ggplot2::aes(x=GLBrH,y=HBRO3))+
    ggplot2::geom_abline(intercept=0,slope=1,col="red")
dd

ddde <- data.frame(pH=pHRES,PKBR=PKBr,BRO3=SPECS[,5],HBRO3=SPECS[,4])

gg <- ggplot2::ggplot()
gg <- gg + ggplot2::geom_line(data=ddde,ggplot2::aes(x=PKBr,y=pH),col="green",
                              size=1)
gg <- gg + ggplot2::geom_line(data=ddde,ggplot2::aes(x=PKBr,y=BRO3),col="red",
                              size=1)
gg <- gg + ggplot2::geom_line(data=ddde,ggplot2::aes(x=PKBr,y=HBRO3),col="blue",
                              size=1)
gg <- gg + ggplot2::geom_point(data=ddde,ggplot2::aes(x=PKBr,y=pH),col="green",
                               size=2)
gg <- gg + ggplot2::geom_point(data=ddde,ggplot2::aes(x=PKBr,y=BRO3),col="red",
                               size=2)
gg <- gg + ggplot2::geom_point(data=ddde,ggplot2::aes(x=PKBr,y=HBRO3),
                               col="blue",size=2)+
    ggplot2::scale_y_continuous(breaks=seq(-0.02,0.12,by=0.02))+
    ggplot2::scale_x_continuous(breaks=seq(-2.5,2.5,by=0.5))+
    ggplot2::xlab("pKBr")+
    ggplot2::ylab("pH (green) and concentrations M")

gg
#Charge balances with Glaser data using their concentrations
H <- brom$H
CHH <- Na + H - kw/H - brom$BrO3-FSPECS[,2]-2*FSPECS[,3]
summary(CHH)

#Charge balances with Glaser data with our methods
H <- 10^-UNCOR_pH
CHH1 <- Na + H - kw/H - FSPECS[,5]-FSPECS[,2]-2*FSPECS[,3]
summary(CHH1)

}
\keyword{datasets}
