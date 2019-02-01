\name{Tessman}
\docType{data}
\alias{Tessman}
\title{Data from Tessman and Ivanov 2002}
\description{Data from table 1 in the paper including mixture of
    acetate, phosphoric acid, boric acid and NaOH}
\usage{Tessman}
\format{WA is list of the buffers with KA values, Na and TOTAL: buffer
    concentrations}
\source{Tessman AB, Ivanov AV. J Analytical chemistry 2002; 57: 2-7}
\examples{data(Tessman)
TOT <- TOTAL
TOTAL <- matrix(TOTAL[1,],c(1,3))
K<-Cl<-Lact<-Ca<-Mg<-Alb<-PCO2<-0
A <- 0.5108;b <- 0.1
Na1 <- Na[1]
pH <- pH_general(Na=Na1,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$pH
IS <- pH_general(Na=Na1,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$II
#print(paste("First point, pH",pH, "\n")
#print(paste("First point, Ionic strength",II, "\n")

#pH
pH
#Ionic strengthh
IS

TOTAL <- TOT
TOTAL <- matrix(TOTAL[2,],c(1,3))
K<-Cl<-Lact<-Ca<-Mg<-Alb<-PCO2<-0
A <- 0.5108;b <- 0.1
Na2 <- Na[2]
(FF<-pH_general(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$FF1)
    print(FF)
    pHA <- pH_general(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$pH
    IS <- pH_general(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$II
    FIS <- pH_general(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$FII
    UNpH <- pH_general(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$UNCOR_pH
    SPECS <- pH_general(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA, TOTAL[1,],
               Alb=Alb,PCO2=PCO2,A=A,b=b)$SPECS
#Second pH
pHA
#Second ionic strength
IS
#Uncorrected pH
UNpH
#First Ionic strength
FIS

    pH <- seq(1,14,length=1000)
    CC <- c()
    for (i in 1:1000) CC[i] <- CMB(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact, WA,
    FF,TOTAL[1,],Alb=Alb,PCO2=PCO2,H=10^-pH[i])

    ddf <- data.frame(pH,CC)
    gg <- ggplot2::ggplot()
    gg <- gg + ggplot2::geom_line(data=ddf,ggplot2::aes(x=pH,y=CC),col="red",
                                  size=1.5)+
        ggplot2::scale_x_continuous(breaks=seq(1,14,by=1)) +
        ggplot2::ylab("Charge M")+
        ggplot2::geom_hline(yintercept=0)+ggplot2::geom_vline(xintercept=pHA)
    print(gg)

    pHS <- round(pHA-0.2,digits=1)
    pHE <- round(pHA+0.2,digits=1)
    print(paste("start pH ", pHS))

    pH1 <- seq(pHS,pHE,length=1000)
    CC1 <- c()

    for (i in 1:1000) CC1[i] <- CMB(Na=Na2,K=K,Cl=Cl,Ca=Ca,Mg=Mg,Lact=Lact,
                                    WA,FF,
                                   TOTAL[1,],Alb=Alb,PCO2=PCO2,H=10^-pH1[i])

    ddf1 <- data.frame(pH1,CC1)
    str(ddf1)
    gg1 <- ggplot2::ggplot()
    gg1 <- gg1 + ggplot2::geom_line(data=ddf1,ggplot2::aes(x=pH1,y=CC1),
                                    col="red",
                                    size=1.5)+
        ggplot2::scale_x_continuous(breaks=seq(pHS,pHE,by=0.05)) +
        ggplot2::ylab("Charge M")+ggplot2::xlab("pH")+
        ggplot2::geom_hline(yintercept=0)+ggplot2::geom_vline(xintercept=pHA)
    print(gg1)

    #Charge balance
    ZZ <- c(0,-1,0,-1,-2,-3,0,-1,-2,-3)

CCH <-  sum(SPECS*ZZ)+0.1+10^-pHA-kw/(FF[1]^2*10^-pHA)
#Our charge balance
CCH
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxx
#  find charge balance at Tessmans equlibrium

Pk1 <- WA$KA[[2]][1];Pk2 <- WA$KA[[2]][2];Pk3 <- WA$KA[[2]][3]
bk1 <- WA$KA[[3]][1];bk2 <- WA$KA[[3]][2];bk3 <- WA$KA[[3]][3]
ak1 <- WA$KA[[1]]

A <- 0.5108; b <- 0.1
FFG <- getfs(c(1,2,3),0.1367,A,b)
TTC <- -(sum(Bspecif(0.02,12.005,c(Pk1,Pk2,Pk3),FFG)[2:4]*c(1,2,3))+
Bspecif(0.02,12.005,ak1,FFG)[2]+
sum(Bspecif(0.02,12.005,c(bk1,bk2,bk3),FFG)[2:4]*c(1,2,3)))+0.1+10^-12.005 -
kw/(FFG[1]^2*10^-12.005) # -0.02235989

TTC/CCH #-330295.5

#total concentrations
(ACETAT <- sum(SPECS[1:2])) #0.02
(PHOSP  <- sum(SPECS[3:6])) #0.02
(BORIC  <- sum(SPECS[7:10])) #0.02

}
\keyword{datasets}

