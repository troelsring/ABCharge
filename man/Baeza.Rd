\name{Baeza}
\alias{Baeza}
\docType{data}
\title{
Data from Baeza-Baeza 2012
}
\description{
Baeza data J Chem Educ 2012; 89: 900-904
}
\usage{data("Baeza")}
\format{
  A workspace.
  \describe{
    \item{Na,K,Cl,Ca,Mg,Lact}{Strong ions.}
    \item{Alb,PCO2}{Albumin (g/l) and pCO2 (mmHg).}
    \item{WA}{ List of buffer names and pks.}
    \item{TOTAL}{Matrix of buffer concentrations.}
  }
}
\details{
 WA TOTAL and Na to reproduce figure 1 in Baeza-Baeza
}
\source{
    Baeza-Baeza JJ et al. J Chem educ 2912; 89: 900-904
}

\examples{
data(Baeza)
j <- 1
pHRES <- c();KpH <- c(); Charge <- IonStr <- iterations <- CCCH <- RCCH <-
UNCOR_pH <- III <- c()
FF1 <- rep(NA,2)
SPECS <- rep(NA,3)
pHRES[j]  <-  pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$pH
UNCOR_pH[j] <- pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$UNCOR_pH
SPECS <- pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$SPECS
FF1  <- pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$FF1
Charge[j]<- pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$Charge
III[j]  <- pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$II
(KpH[j]  <- pH_general(Na=Na,K=K,Cl=Cl,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=WA,
TOT=TOTAL[j,],A=0.509,b=0.1)$KpH)


cat("Our corected pH")
KpH

#Charge balance from Baeza figure 1
CHHBA <- 9.796e-4+5e-2-1.88e-2-3.16e-2-2*2.911e-4-1.555e-11 #-2.600016e-06
#Phtalate  total concentratrion
PHT <- 3.16e-2 + 1.81e-2 + 2.911e-4 #0.0499911
#Our charge balance
Charge
#Our total Pth
sum(SPECS[1:3])

#Recovery of pKs

(RPK1 <- -log10(FF1[1]^2*10^-pHRES*SPECS[2]/SPECS[1])) #2.95000010687057
(RPK2 <- -log10(FF1[2]*10^-pHRES*SPECS[3]/SPECS[2]))   #5.41000248868762

#Baeza
(R1B <- -log10((10^-0.091)^2*9.796e-4*3.160e-2/1.811e-2)) #2.94918259124298
(R2B <- -log10((10^-0.365)*9.796e-4*2.911e-4/3.160e-2))   #5.40959610072692
}
\keyword{datasets}
