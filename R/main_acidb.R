
#' Find equilibrium parameters
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#  THIS IS THE MAIN ENTRY
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#' @export
pH_general <- function(Na=0,K=0,Cl=0,Lact=0,Ca=0,Mg=0,PCO2=0,Alb=0,WA=list(),TOT,A=A,b=b) {


#' @param Na,K,Cl,Lact,Ca,Mg  Strong ions inputs, "0" if not present.
#' @param PCO2,Alb If present PCO2 and Albumin are handled - "0" if not present.
#' @param WA List of lists with names of all buffers in first list, and all pK
#'  values from smallest to highest for each buffer in second list.
#' @param TOT Vector of total buffer concentrations.
#' @param A,b Constants from Davies equation.
#' @return Returns pH: pH as -log10 of direct H from charge balance
#'   minimization after correction for ionic strength, KpH: -log10 of that H
#'   multiplied by F1 (activity coeff for charge=1), UNCOR_pH: pH with all
#'   activity coeffs = 1,II: final ionic strength, FII: ionic strength before
#'   optimizing activity coefficients, i: number of iterations in optimization,
#'   SPECS: vector of all buffer species, FSPECS: species distribution before
#'   optimizing activity coefficients.
#' @importFrom stats uniroot


  absdiff <- 1e-5 #pH difference between iterations
  diff <- 1
  i <- 1  #iteration number

  FF <- rep(1,6)


  while(diff>absdiff) {


    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #                FIND first pH
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx





    pH_START <-  -log10(uniroot(CMB,c(1e-20,5),tol=.Machine$double.eps,maxiter=100000,
                                Na=Na,K=K,Ca=Ca,Mg=Mg,Cl=Cl,Lact=Lact,FF=FF,TOT=TOT,WA=WA,Alb=Alb,PCO2=PCO2)$root)


    if(i==1) {UNCOR_pH <- pH_START}

    Charge <- CMB(Na=Na,K=K,Ca=Ca,Mg=Mg,Cl=Cl,Lact=Lact,FF=FF,TOT=TOT,WA=WA,Alb=Alb,PCO2=PCO2,H=10^-pH_START)


    #GET_CH(1,FF,WA,10^-4.28)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #  Based on WA and pH and TOT find specific moities via Bspecif(TTOT,pH,KA) - gather them in vector BUFFS
    # here again is an issue if KA is only 1 item
    # Also here: all species recovered
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    BUFFS <-  LL <- SPECS <- c()
    for (N in 1:length(WA$buffs))
    {
      BUFFS <- c(BUFFS,Bspecif(TOT[N],pH_START,WA$KA[[N]],FF)[2:(length(WA$KA[[N]])+1)])
      SPECS <- c(SPECS,Bspecif(TOT[N],pH_START,WA$KA[[N]],FF))
      LL[N]  <-length(WA$KA[[N]])
    }


    if (i==1) FSPECS = SPECS
    CCC <- c() #Find charges
    for (j in 1:length(LL)) CCC <- c(CCC,seq(-1,-LL[j]))
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Include CO2 species in ionic strength assessment
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    CCO2 <- c()
    CO2Z <- c()
    if(PCO2>0) {CCO2 <- CO2specif(PCO2,pH_START,FF) #returns c(HCO3,CO3) if PCO2 is present
    CO2Z <- c(-1,-2)}        #returns relevant charges if needed


    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Make CONCentration and charge vectors for calculating ionic strength - exclude albumin
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    CONC <- c(Na,Cl,K,Ca,Mg,Lact,BUFFS,10^-pH_START,kw/(10^-pH_START),CCO2)
    ZZZ <- c(1,-1,1,2,2,-1,CCC,1,-1,CO2Z)
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    II <- getI(CONC,ZZZ)
    # print(II)
    if (i == 1) FII <- II


    FFN <-getfs(seq(1,max(abs(ZZZ))),II,A,b)
    FF <- FFN   #FF updated


    pH_NEW <-  -log10(uniroot(CMB,c(1e-20,5),tol=.Machine$double.eps,maxiter=100000,
                              Na=Na,K=K,Ca=Ca,Mg=Mg,Cl=Cl,Lact=Lact,FF=FF,TOT=TOT,WA=WA,Alb=Alb,PCO2=PCO2)$root)


    diff  <- abs(pH_NEW-pH_START)

    i <- i +1

    KpH <- -log10(FF[1]*10^-pH_NEW)

  }
  return(list(pH= pH_NEW,KpH=KpH,i=i,diff=diff,II = II,Charge=Charge,FF1=FF,UNCOR_pH=UNCOR_pH,FII=FII,SPECS=SPECS,FSPECS=FSPECS))
}

