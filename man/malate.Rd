\name{malate}
\alias{malate}
\alias{pHOBS}
\docType{data}
\title{
Malate
}
\description{
Malic acid data from Glaser 2014
}
\usage{data("malate")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{pHOBS}{Observed pH}
     }
}
\details{
Malic acid and malate data fro Glaser 2014
}
\source{
Glaser et al Dynamical approach  J Chem Eduuc 2014; 91: 1009-1016
}
\examples{
data(malate)
pHRES <- c();KpH <- c(); Charge <- IonStr <- iterations <- CCCH <-
RCCH <- UNCOR_pH <- III <- c()
FF1 <- array(rep(NA,2*33),c(33,2))
SPECS <- array(rep(NA,3*33),c(33,3))

for (j in 1:33) {
  pHRES[j]  <-  pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$pH
  UNCOR_pH[j] <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$UNCOR_pH
  SPECS[j,] <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$SPECS
  FF1[j,]  <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$FF1
  Charge[j]<- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$Charge
  III[j]  <- pH_general(Na=Na[j],K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,
  WA=WA,TOT=TOTAL[j,],A=0.5108,b=0.1)$II
}

ddf <- data.frame(pHC=pHRES,pHO=UNCOR_pH,II=III,pHOBS=pHOBS)

ss <- ggplot2::ggplot()
ss <- ss + ggplot2::geom_point(data=ddf,ggplot2::aes(x=pHOBS,y=pHO))
ss <- ss + ggplot2::geom_abline(intercept=0,slope=1,col="red")+
ggplot2::xlab("Observed pH") + ggplot2::ylab("Modeled pH")
ss + ggplot2::ggtitle("Malate reults from Glaser 2014
             observed versus modeled pH assuming ideal state")

#Correlation
(XX <- cor.test(pHOBS,UNCOR_pH)$estimate)

gg <- ggplot2::ggplot()
gg <- gg + ggplot2::geom_point(data=ddf,ggplot2::aes(x=pHO,y=pHC,col=II))+
ggplot2::xlab("pH ignoring ionic strength")+ggplot2::ylab(
"pH corrected for ionic strength")+
  ggplot2::geom_abline(intercept=0,slope=1,col="red")+
  ggplot2::scale_colour_continuous(name="Ionic strength")

gg

dd <- ggplot2::ggplot()
dd <- dd + ggplot2::geom_point(data=ddf,ggplot2::aes(x=II,y=pHO-pHC))
dd + ggplot2::xlab("Ionic strength") + ggplot2::ylab("pH difference") +
  ggplot2::ggtitle("Difference between pH assuming ideality and corrected pH
          as function of ionic strength")

}
\keyword{datasets}
