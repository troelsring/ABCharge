

#' @export
kw <- 1e-14 # 2.39e-14
#' @export
Kc <- 2.45e-11
#' @export
K3 <- 5.76e-11
KH <- 1.77e-7
AF <- 21/66500 # fixed negatives
AH <- 16/66500 # Histidines


#' BCharge

#' Finds charge on any buffer given activity coefficients (FF), H, KA,tot

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# Finds charge on any buffer given activity coefficients (FF), H, KA,tot
# #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @export
BCharge <- function(TOT,FF,KA,H)
{
  #' @param TOT vectore of buffer concentrations
  #' @param KA vector of dissociation constants
  #' @param FF vactor of activity coefficients
  #' @param H  proton concentration
  #'
  num <- c()
  num[1] <- KA[1]/(FF[1]^2*H)
  if (length(KA)>1){
    for (i in 2:length(KA)) num[i] <- i*prod(KA[1:i])/(FF[1]^i*FF[i]*H^i)}
  num <- sum(num)*TOT
  denum <- c()
  denum[1] <- 1+ KA[1]/(FF[1]^2*H)
  if (length(KA)>1){
    for (i in 2:length(KA)) denum[i] <- prod(KA[1:i])/(FF[1]^i*FF[i]*H^i)}
  denum <- sum(denum)
  num/denum
}

#' GET_CH

#' Calls BCharge to gather charges from a buffer
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#Sums the buffer charges at given TOT,FF,H and at given  KA's
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#' @export
GET_CH <- function(TOT,FF,WA,H){
    #' @param TOT vectore of buffer concentrations
    #' @param WA list of buffer specifications
    #' @param FF vactor of activity coefficients
    #' @param H  proton concentration
  CC <- c()
  for (i in 1: length(WA$KA))   CC[i] <- BCharge(TOT[i],FF,KA=WA$KA[[i]],H)
  sum(CC)
}

#' CMB
#' Generic function to find charge

#
# xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#             FINDS CURRENT pH
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#' @export
CMB  <- function(Na,K,Cl,Ca,Mg,Lact,H,WA,FF,TOT,PCO2,Alb)
    #'

    #' @param Na,K,Cl,Ca,Mg,Lact Strong ions for the charge balance
    #' @param WA,TOT  Buffer specifications and total concentrations
    #' @param FF,H,PCO2,Alb  Activity coefficients, H, PCO2 and albumin

    {Na+K+2*Ca+2*Mg-Cl-
    Lact+H-kw/(FF[1]^2*H)-GET_CH(TOT,FF,WA,H)-
    PCO2*(Kc/(FF[1]^2*H)+(2*K3*Kc)/(FF[1]^3*H^2*FF[2]))-AF*Alb+AH*Alb*H/(KH+H)
}


#' Bspecif
#' Finds the partition of buffer species
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#   FINDS THE PARTITION OF BUFFER SPECIES - must be called with each non CO2
# buffer via KA, the vector of KA, and TTOT concentration
#   Here we also need care if KA is only 1 item
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @param TTOT,pH Total buffer concentration and pH
#'
#' @param KA,FF vectors of dissociations and activivty coefficients

#' @export
Bspecif  <- function(TTOT,pH,KA,FF) {

  H <- 10^-pH
  denum <- c()

  denum[1] <- 1+ KA[1]/(FF[1]^2*H)
  if (length(KA)>1){
    for (i in 2:length(KA)) denum[i] <- prod(KA[1:i])/(FF[1]^i*FF[i]*H^i)

  }
  denum <- sum(denum)
  HX <- TTOT/denum  #Uncharged fully H loaded moiety
  num <- c()
  for (i in 1:length(KA))num[i]<- HX*(prod(KA[1:(i)])/(FF[1]^(i)*FF[(i)]*H^(i)))
  num <- c(HX,num)
  num
}

#' CO2specif - finds the partition between CO2 species
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#      FINDS THE PARTITION BETWEEN BICRBONATE AND CARBONATE
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#'
#' @param PCO2,pH,FF  CO2 tension, pH and activity coefficients
#' @export
CO2specif  <- function(PCO2,pH,FF) {
  H <- 10^-pH

  HCO3 <- PCO2*Kc/(FF[1]^2*H)
  CO3  <- PCO2*K3*Kc/(FF[1]^3*H^2*FF[2])
  c(HCO3,CO3)
}


#' GetI - finds the ionic strength
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#     FINDS IONIC STRENGTH GIVEN ALL CONCENTRATIONS FROM INPUT OF STRONG
#    IONS AND ALL BUFFERS
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#' @export
getI <- function(CONC,VALS)
    #' @param CONC,VALS  concentrations and valencies
{0.5*sum(VALS^2 * CONC)}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#            GET ALL ACTIVITY COEFFS CORRESPONDING TO z
# GIVEN I AND b, USUALLY 0.1
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Getfs
#'
#' Gets acitivity coefficiens from Davies equation.

#' @export
getfs <- function(z,I,A,b) {
    #' @param A,b are parameters in Davies equation
    #' @param z,I - vector of charges and ionic strength

  10^(-A*z^2*(sqrt(I)/(1+sqrt(I))-b*I))}

