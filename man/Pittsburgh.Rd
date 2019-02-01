\name{Pittsburgh}
\alias{Pittsburgh}
\alias{INPUT}
\alias{WA}
\alias{PCO2}
\alias{Alb}
\alias{Na}
\alias{K}
\alias{Cl}
\alias{Ca}
\alias{Mg}
\alias{Lact}

\docType{data}
\title{
Pittsburgh ICU data
}
\description{
2437 measurements from ICU
}
\usage{data("Pittsburgh")}
\format{
  A workspace with clinical data.
  \describe{
    \item{Na}{Strong ion.}
    \item{K}{Strong ion.}
    \item{PCO2}{Tension of pCO2 in mmHg.}
    \item{INPUT}{Dataframe of data from ICU.}
    \item{INPUT$Cl}{Cl in INPUT}
  }
}
\source{
Ring T, Kellum JA. Strong relationships in acid-base chemistry – modeling protons based on predictable concentra-tions of strong ions, total weak acid concentrations, and pCO2. Plos One 2016 Sep 15;11(9):e0162872. doi: 10.1371/journal.pone.0162872. eCollection 2016.PMID: 27631369
}
\examples{
data(Pittsburgh)

pH_obs <- INPUT$pH



pHRES <- c();KpH <- c(); IonStr <- iterations <- CCCH <- RCCH <- II <-
UNCOR_pH <- c()
FFF <- array(rep(NA,6*length(INPUT[,1])),c(length(INPUT[,1]),6))

for (j in 1:length(INPUT[,1]))

  {

  if (NCOL(TOTAL)==1) {
    TOT <- TOTAL[j]
  } else {TOT <-TOTAL[j,] }
  pHRES[j] <- with(INPUT,pH_general(Na=Na[j]*1e-3,K=K[j]*1e-3,Cl=Cl[j]*1e-3,
  Lact=Lact[j]*1e-3,Ca=Ca[j],Mg=Mg[j],PCO2=PCO2[j],
                                    Alb=ALB[j],WA=WA,TOT=TOT,A=0.5108,b=0.1)$pH)
  FFF[j,] <- with(INPUT,pH_general(Na=Na[j]*1e-3,K=K[j]*1e-3,Cl=Cl[j]*1e-3,
  Lact=Lact[j]*1e-3,Ca=Ca[j],Mg=Mg[j],PCO2=PCO2[j],
                                   Alb=ALB[j],WA=WA,TOT=TOT,A=0.5108,b=0.1)$FF1)
  KpH[j] <- with(INPUT,pH_general(Na=Na[j]*1e-3,K=K[j]*1e-3,Cl=Cl[j]*1e-3,
  Lact=Lact[j]*1e-3,Ca=Ca[j],Mg=Mg[j],PCO2=PCO2[j],
                                  Alb=ALB[j],WA=WA,TOT=TOT,A=0.5108,b=0.1)$KpH)
   II[j] <- with(INPUT,pH_general(Na=Na[j]*1e-3,K=K[j]*1e-3,Cl=Cl[j]*1e-3,
   Lact=Lact[j]*1e-3,Ca=Ca[j],Mg=Mg[j],PCO2=PCO2[j],
                                  Alb=ALB[j],WA=WA,TOT=TOT,A=0.5108,b=0.1)$II)
   UNCOR_pH[j] <- with(INPUT,pH_general(Na=Na[j]*1e-3,K=K[j]*1e-3,Cl=Cl[j]*1e-3,
   Lact=Lact[j]*1e-3,Ca=Ca[j],Mg=Mg[j],PCO2=PCO2[j],
                            Alb=ALB[j],WA=WA,TOT=TOT,A=0.5108,b=0.1)$UNCOR_pH)
}
#ionic strength
summary(II)

summary(z5 <- lm(KpH~pH_obs))

dd <- data.frame(pH_obs,KpH)

XXX <- cor.test(pH_obs,KpH)$estimate

gg <- ggplot2::ggplot()
gg <- gg + ggplot2::geom_point(data=INPUT,ggplot2::aes(x=pH_obs,y=KpH))+
ggplot2::xlab("Measured pH")+ggplot2::ylab("Modeled pH")+
    ggplot2::geom_abline(intercept=0,slope=1,col="red",size=2)+
  ggplot2::geom_abline(intercept=0.69415,slope=0.90319,col="blue",size=2)+
  ggplot2::ggtitle("Red line is x=y, blue line is fitted model\n
          modeled pH is based on Davies activity coefficient")+
  ggplot2::scale_x_continuous(limits=c(6.8,7.8))+
  ggplot2::scale_y_continuous(limits=c(6.8,7.8))
gg
ddif1 <- KpH-pH_obs
ddif2 <- pHRES-pH_obs
ddif3 <- UNCOR_pH-pH_obs
summary(ddif1)
summary(ddif2)
summary(ddif3)

grp <-  c(rep("blue",length(ddif1)),rep("red",length(ddif1)),
rep("green",length(ddif1)))
ssd <- data.frame(ddif=c(ddif1,ddif2,ddif3),col=grp)

ggplot2::ggplot(data=ssd,ggplot2::aes(x=ddif,fill=col))+
ggplot2::geom_density(alpha=0.2)+ggplot2::xlab("Difference: modeled-observed")+
  ggplot2::labs(fill="Group")+
  ggplot2::scale_fill_discrete(labels=c("Corrected pH difference",
  "pH difference assuming ideal conditions","Uncorrected pH difference"))+
  ggplot2::geom_vline(xintercept=0,linetype=2,size=1)

}

\keyword{datasets}
