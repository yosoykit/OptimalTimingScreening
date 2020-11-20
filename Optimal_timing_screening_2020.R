#######################################################
## OPTIMAL Timing PAPER RESULTS 2020 November        ##
##      KIT CURTIUS ET AL.                           ##
##   Computes Figure 2, Table 1, Costs,              ##
##  Figures 3-5, and Supplementary Fig 1B.           ##
##   Some scripts in other Github repos as indicated ##
#######################################################

source('/BEtoEAC_Results/am_parameter_list_2020.R')
source('BE_density_Optimal_R.R')
source('BaseCaseII/AllMortalityGenUS.R')
source('/BaseCaseII/legauss.R')
dyn.load('/BaseCaseII/legauss.so')
library('deSolve')
library('fields')

##################################################################################
##### FUNCTIONS FOR ANALYTICAL SOLUTIONS TO TIMING EQUATIONS GIVEN BELOW #########
##################################################################################
phidot = function(s, phi, parms){
  phidot=numeric(10)
  with(as.list(parms),{
  tau = t - s
  #RR=5
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))
  nu = nu0*((1-GERDpr)+RR*GERDpr)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## A-D transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}


phidot_GERD = function(s, phi, parms){
  phidot=numeric(10)
  with(as.list(parms),{
  tau = t - s
  #RR=5
  GERDpr=1-exp(-gerd1*min(gerd3,tau)-gerd2*max(0,tau-gerd3))
  #nu = nu0*((1-GERDpr)+RR*GERDpr)
  nu = nu0*(RR*1)
  phidot[1]=betam -(alpham+betam+rho)*phi[1]+alpham*(phi[1]**2)
  phidot[2]=2*alpham*phi[1]*phi[2]- (alpham+betam+rho)*phi[2]
  phidot[3]=betap+ mu2*phi[1]*phi[3]-(alphap+betap+mu2)*phi[3]+alphap*(phi[3]**2)
  phidot[4]=2*alphap*phi[3]*phi[4]+mu2*(phi[4]*phi[1]+phi[3]*phi[2]) - (alphap+betap+mu2)*phi[4]
  phidot[5]=mu1*phi[5]*(phi[3]-1)
  phidot[6]=mu1*(phi[6]*(phi[3]-1)+phi[5]*phi[4])
  phidot[7]=mu0*phi[7]*(phi[5]-1)
  phidot[8]=mu0*(phi[8]*(phi[5]-1)+phi[7]*phi[6])
  ## A-D transition
  phidot[9]=nu*(phi[7]-phi[9])
  phidot[10]=nu*(phi[8]-phi[10])
  list(c(phidot))
  })
}

s_EAC <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi9[2])
   }
  return(besurv.ode)
}  

h_EAC <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,-ode$phi10[2]/ode$phi9[2])
   }
  return(besurv.ode)
}  

s_EAC_GERD <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi9[2])
   }
  return(besurv.ode)
}  

h_EAC_GERD <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,-ode$phi10[2]/ode$phi9[2])
   }
  return(besurv.ode)
}  

s_MSCE <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi7[2])
   }
  return(besurv.ode)
}  
s_3 <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi5[2])
   }
  return(besurv.ode)
}  
s_2 <- function(parms,t){
  ages<- t
  besurv.ode = NULL
  for (a in ages) {
      # initial condition
      parms['t'] = a
      phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0,phi9=1,phi10=0)
      times = c(0,a)
      ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
      besurv.ode = c(besurv.ode,ode$phi3[2])
   }
  return(besurv.ode)
}  

p_2 <- function(parms,ts){
  with(as.list(parms),{
  zeta = (exp((alphap-betap)*t)-1)/(alphap*(exp((alphap-betap)*t)-betap))
  return(1-betap*zeta)
  })
}


###################################
###     SCREENING  Timing       ###
###################################


#########################################################
### Screening strategy 1: RESULTS FOR MALES BORN 1950 ###
#########################################################
RR=5
b_c1 = 1950
start=10
w=seq(0,10, .25)
t_s=start:80
birth_cohort = b_c1
survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv[1:length(t_s)]
survdeath_by_OCf <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv[1:length(t_s)]

source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]


p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))
strat1_w_GERD = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

for (i in 1:length(w)){
  strat1_w[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OC
  strat1_w_GERD[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OC
}



###########################################################
### Screening strategy 1: RESULTS FOR FEMALES BORN 1950 ###
###########################################################

source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)
nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,nu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,nu0,RR)$nu_cum[t_s]

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))
strat1_w_GERD_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

for (i in 1:length(w)){
  strat1_w_f[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OCf
  strat1_w_GERD_f[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OCf
}


###################################
###     Figure 2 (A-D)          ###
###################################
dev.new()
par(mfrow=c(2,2), mai = c(0.2, 0.05, 0.02, 0.1))
z=strat1_w
y=w
x=t_s
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("purple", "blue", 'green') )
# Generate the desired number of colors from this palette
nbcol <- 10
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
col_range=seq(range(z)[1],range(z)[2],(range(z)[2]-range(z)[1])/nbcol)

res=persp(x, y, z, col = color[facetcol], phi = 30, theta = 30,xlab="Screening time (age)", ylab="\n\n\n w (weight parameter)",expand=.8,ticktype = "detailed",zlab="\n\n\n\n Screening objective function",border=NA, cex.axis=1.25, cex.lab=1.25)
image.plot(legend.only=T, zlim=range(z), col=color,add=T,legend.mar=15)
clines=contourLines(x,y,z)
lapply(clines,function(contour){lines(trans3d(contour$x, contour$y,range(z)[1],res),col=color[which.min(abs(col_range - contour$level))])})
lines(trans3d(t_s[apply(strat1_w,2,which.max)],w,z=range(z)[1], pm=res),col='red',lwd=1,type='l', cex=.5)
points(trans3d(t_s[apply(strat1_w,2,which.max)][5],1,z=range(z)[1], pm=res),col='red', cex=1.5,pch=17)
#title(main="A)     All Males")

#########################################################
###  RESULTS FOR GERD MALES BORN 1950                 ###
#########################################################

z=strat1_w_GERD
y=w
x=t_s
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("purple", "blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 10
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
col_range=seq(range(z)[1],range(z)[2],(range(z)[2]-range(z)[1])/nbcol)
res=persp(x, y, z, col = color[facetcol], phi = 30, theta = 30,xlab="Screening time (age)", ylab="\n\n\n w (weight parameter)",expand=.8,ticktype = "detailed",zlab="\n\n\n\n Screening objective function",border=NA, cex.axis=1.25, cex.lab=1.25)
image.plot(legend.only=T, zlim=range(z), col=color,add=T,legend.mar=15)
clines=contourLines(x,y,z)
lapply(clines,function(contour){lines(trans3d(contour$x, contour$y,range(z)[1],res),col=color[which.min(abs(col_range - contour$level))])})
lines(trans3d(t_s[apply(strat1_w_GERD,2,which.max)],w,z=range(z)[1], pm=res),col='red',lwd=1,type='l', cex=.5)
points(trans3d(t_s[apply(strat1_w_GERD,2,which.max)][5],1,z=range(z)[1], pm=res),col='red', cex=1.5,pch=17)
#title(main="B)     Males with GERD symptoms")

z=strat1_w_f
y=w
x=t_s
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("purple", "blue", 'green') )
# Generate the desired number of colors from this palette
nbcol <- 10
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
col_range=seq(range(z)[1],range(z)[2],(range(z)[2]-range(z)[1])/nbcol)

res=persp(x, y, z, col = color[facetcol], phi = 30, theta = 30,xlab="Screening time (age)", ylab="\n\n\n w (weight parameter)",expand=.8,ticktype = "detailed",zlab="\n\n\n\n Screening objective function",border=NA, cex.axis=1.25, cex.lab=1.25)
image.plot(legend.only=T, zlim=range(z), col=color,add=T,legend.mar = 15)
clines=contourLines(x,y,z)
lapply(clines,function(contour){lines(trans3d(contour$x, contour$y,range(z)[1],res),col=color[which.min(abs(col_range - contour$level))])})
lines(trans3d(t_s[apply(strat1_w_f,2,which.max)],w,z=range(z)[1], pm=res),col='red',lwd=1,type='l', cex=.5)
points(trans3d(t_s[apply(strat1_w_f,2,which.max)][5],1,z=range(z)[1], pm=res),col='red', cex=1.5,pch=19)
#title(main="C)     All Females")

#########################################################
###  RESULTS FOR GERD FEMALES BORN 1950               ###
#########################################################

z=strat1_w_GERD_f
y=w
x=t_s
nrz <- nrow(z)
ncz <- ncol(z)
# Create a function interpolating colors in the range of specified colors
jet.colors <- colorRampPalette( c("purple", "blue", "green") )
# Generate the desired number of colors from this palette
nbcol <- 10
color <- jet.colors(nbcol)
# Compute the z-value at the facet centres
zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
# Recode facet z-values into color indices
facetcol <- cut(zfacet, nbcol)
col_range=seq(range(z)[1],range(z)[2],(range(z)[2]-range(z)[1])/nbcol)
res=persp(x, y, z, col = color[facetcol], phi = 30, theta = 30,xlab="Screening time (age)", ylab="\n\n\n w (weight parameter)",expand=.8,ticktype = "detailed",zlab="\n\n\n\n Screening objective function",border=NA, cex.axis=1.25, cex.lab=1.25)
#lines(trans3d(t_s[apply(strat1_w,2,which.max)],w,z=apply(strat1_w,2,max), pm=res),col='red',lwd=1,type='o',pch=8, cex=.5)
image.plot(legend.only=T, zlim=range(z), col=color,add=T,legend.mar=15 )
clines=contourLines(x,y,z)
lapply(clines,function(contour){lines(trans3d(contour$x, contour$y,range(z)[1],res),col=color[which.min(abs(col_range - contour$level))])})
lines(trans3d(t_s[apply(strat1_w_GERD_f,2,which.max)],w,z=range(z)[1], pm=res),col='red',lwd=1,type='l', cex=1)

points(trans3d(t_s[apply(strat1_w_GERD_f,2,which.max)][5],1,z=range(z)[1], pm=res),col='red', cex=1.5,pch=19)
#title(main="D)     Females with GERD symptoms")






#########################################################################
### MAIN RESULTS FOR VALIDATION WITH CORI FOR MALES/FEMALES BORN 1950 ###
#########################################################################

w=seq(0,10, .25)

source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}



#########################################################
### Screening setup: RESULTS FOR MALES BORN 1950      ###
#########################################################

b_c1 = 1950
start=10
w=seq(0,10, .25)
t_s=start:80
RR=5
birth_cohort = b_c1
survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv[1:length(t_s)]

source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]


p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))
strat1_w_GERD = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

for (i in 1:length(w)){
  strat1_w[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OC
  strat1_w_GERD[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OC
}



######################################################################################################
###       RESULTS FOR SCREENING STRATEGIES AND ASSOCIATED EFFICACY METRICS                         ###
###                             Table 1                                                            ###
######################################################################################################

################################################################
###    COMPUTE Successful Diagnoses Function (SD(t_s^*)):     ###
##   Pr[ T_BE <= t | t_s^* < T_EAC < b],   0 <= t <= t_s^*   ###
################################################################

#################################
### TABLE 1.  FOR MALES.  #####
#################################
RR=5
source('/BEtoEAC_Results/am_parameter_list_2020.R')

b=80

parms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho =allmales$rho,t=0)
parms['t'] = b
phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,b)
ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
surv80 <- ode$phi9[2]

surv=surv80

ts_star=64
parms['t'] = ts_star
phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,ts_star)
ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
surv_ts_star <- ode$phi9[2]
surv_ts_star_allmales = surv_ts_star

final_int=b
gauss_points = 20

BEonsetcum_MSCE80<-rep(1,b)
for (i in 1:(ts_star)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    parms['t'] = ts_star-xg[a]
    times = c(0,(ts_star-xg[a]))
    ode_i = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    GERDpr1=1-exp(-gerd1*min(gerd3,xg[a])-gerd2*max(0,xg[a]-gerd3))
    nu <- amnu0*((1-GERDpr1)+RR*GERDpr1)
    GLquad1<-legauss(0 ,xg[a],20)
    wg2 <- GLquad1$weights
    xg2 <- GLquad1$mesh
    int1 <- rep(0,20)
    for ( j in 1:20){
      GERDpr=1-exp(-gerd1*min(gerd3,xg2[j])-gerd2*max(0,xg2[j]-gerd3))
      nu_temp <-nu0*((1-GERDpr)+RR*GERDpr)
      int1[j] <- nu_temp
    }
    nu_cum<-1-exp(-sum(wg2*int1))
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode_i$phi7[2]-ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcum_MSCE80[i] <- (1/(surv_ts_star-surv))*(f_int)
}


## Compute SD probability 
SD_allmales = BEonsetcum_MSCE80[ts_star]


### FOR MALES WITH GERD SYMPTOMS
parms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho =allmales$rho,t=0)
parms['t'] = b
phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,b)
ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
surv80GERD <- ode$phi9[2]

ts_star=58
#ts_star = 50
parms['t'] = ts_star
phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,ts_star)
ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
surv_ts_star <- ode$phi9[2]
surv_ts_star_GERDmales = surv_ts_star

final_int=b
gauss_points = 20

BEonsetcumGERD_MSCE80<-rep(1,b)
for (i in 1:(ts_star)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    parms['t'] = ts_star-xg[a]
    times = c(0,(ts_star-xg[a]))
    ode_i = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    nu <- nu0*(RR*1)
    nu_cum<-1-exp(-nu*xg[a])
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode_i$phi7[2]-ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcumGERD_MSCE80[i] <- (1/(surv_ts_star-surv80GERD))*(f_int)
}

## Compute SD probability 

SD_GERDmales=BEonsetcumGERD_MSCE80[ts_star]

################################################################
###    COMPUTE Overdiagnosis Function (OD(t_s^*)):           ###
##            Pr [ T_BE <= t_s | T_EAC > B]                  ###
################################################################

gauss_points = 20

BEonsetcum_MSCE_noEAC80<-rep(0,b)
for (i in 1:(final_int)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  f4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    GERDpr1=1-exp(-gerd1*min(gerd3,xg[a])-gerd2*max(0,xg[a]-gerd3))
    nu <- amnu0*((1-GERDpr1)+RR*GERDpr1)
    GLquad1<-legauss(0 ,xg[a],20)
    wg2 <- GLquad1$weights
    xg2 <- GLquad1$mesh
    int1 <- rep(0,20)
    for ( j in 1:20){
      GERDpr=1-exp(-gerd1*min(gerd3,xg2[j])-gerd2*max(0,xg2[j]-gerd3))
      nu_temp <-amnu0*((1-GERDpr)+RR*GERDpr)
      int1[j] <- nu_temp
    }
    nu_cum<-1-exp(-sum(wg2*int1))
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcum_MSCE_noEAC80[i] <- (1/(surv80))*(f_int)
}

## Compute Over-screened probability 
ts_star = 64
prob_overdiag_noEAC80 = survdeath_by_OC[(ts_star-9)]*(BEonsetcum_MSCE_noEAC80[ts_star]/strat1_w[(ts_star-9),5])*surv80
OD_allmales = BEonsetcum_MSCE_noEAC80[ts_star]


### GERD males only

BEonsetcumGERD_MSCE_noEAC80<-rep(0,b)
for (i in 1:(final_int)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    nu <- nu0*(RR*1)
    nu_cum<-1-exp(-nu*xg[a])
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcumGERD_MSCE_noEAC80[i] <- (1/(surv80GERD))*(f_int)
}
########################################
## Compute Over-screened probability  ##
########################################

ts_star = 58
prob_overdiagGERD_noEAC80 = survdeath_by_OC[(ts_star-9)]*(BEonsetcumGERD_MSCE_noEAC80[ts_star]/strat1_w_GERD[(ts_star-9),5])*surv80GERD
OD_GERDmales = BEonsetcumGERD_MSCE_noEAC80[ts_star]


##################
## Compute PPV: ##
##################
osGERD = prob_overdiagGERD_noEAC80
ssGERD = 1 - osGERD
os = prob_overdiag_noEAC80 
ss = 1- os
ppv_male = ss
ppv_GERDmale = ssGERD

print(paste("PPV GERD+ males:",ppv_GERDmale))
print(paste("PPV all males:",ppv_male))
print(paste("SD all males:",SD_allmales))
print(paste("SD GERD+ males:",SD_GERDmales))
print(paste("OD all males:",OD_allmales))
print(paste("OD GERD+ males:",OD_GERDmales))
print(paste("Cancer risk all males:",1-surv_ts_star_allmales))
print(paste("Cancer risk GERD+ males:",1-surv_ts_star_GERDmales))


#####################################################
## Compute COST comparison example for GERD+ males ##
#####################################################
################################################
## GERD+ males t_s^*=58  ----> optimal age    ##
########      vs.                             ##
## GERD+ males t_s=50  ----> current practice ##
################################################
source('/BEtoEAC_Results/am_parameter_list_2020.R')

ts_star50=50
parms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho =allmales$rho,t=0)
parms['t'] = ts_star50
phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,ts_star50)
ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
surv_ts_50GERDmales <- ode$phi9[2]

BEonsetcumGERD_MSCE_noEAC80<-rep(0,b)
for (i in 1:(final_int)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    nu <- nu0*(RR*1)
    nu_cum<-1-exp(-nu*xg[a])
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcumGERD_MSCE_noEAC80[i] <- (1/(surv80GERD))*(f_int)
}

osGERD50= survdeath_by_OC[(ts_star50-9)]*(BEonsetcumGERD_MSCE_noEAC80[ts_star50]/strat1_w_GERD[(ts_star50-9),5])*surv80GERD


### COST FUNCTION FOR # of endoscopies (C1) and corresponding csts (C2)

## Number endos / BE patient-year = Number of surveillance endoscopies with RFA HGD strategy / (1000 BE patients * average lifespan in simulation) -> ref: Kroep*, Heberle*, Curtius*, et al. Clin Gastro Hepatol  2017 [18] - FHCRC model Supplementary Table E2
## Number endos / BE patient-year * OS (t_s)*(80 - t_s)*BE_prev(t_s)* Number of initial screens
# t_s = 58
# BE prevalence in cancer-free population undergoing screening: c(t_s[seq(15,65,1)][35],strat1_w_GERD_noMortality50[seq(15,65,1),5][35])
#[1] 58.00000000  0.09706724
ts_star=58
cost1_58=osGERD*(b-ts_star) *strat1_w_GERD_noMortality50[seq(15,65,1),5][35]
print(cost1_58)
#[1] 1.524894
cost2_58 = (7798 /(1000*19.5))*cost1_58
print(cost2_58)
#[1] 0.6098012

# t_s = 50 
# BE prevalence in cancer-free population undergoing screening: c(t_s[seq(15,65,1)][27],strat1_w_GERD_noMortality50[seq(15,65,1),5][27])
# [1] 50.00000000  0.08640053
ts_star50=50
cost1_50= osGERD50*(b-ts_star50)*strat1_w_GERD_noMortality50[seq(15,65,1),5][27]
print(cost1_50)
#[1] 1.67641
cost2_50=(7798 /(1000*19.5))*cost1_50
print(cost2_50)
#[1] 0.6703919

# Number of additional EGDs in overdiagnosed population
print(0.67*400000-.61*400000)
# [1] 24000
# Costs in USD for additional EGDs (assuming EGS costs $745 each. ref: Heberle et al. 2017 [19])
print(24000*745)
#[1] 17880000


#################################
### TABLE 1.  FOR FEMALES.  #####
#################################
source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
  strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
  strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}

###########################################################
### Screening setup: RESULTS FOR FEMALES BORN 1950.     ###
###########################################################
b_c1 = 1950
start=10
w=seq(0,10, .25)
t_s=start:80
RR=5
birth_cohort = b_c1
survdeath_by_OCf <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv[1:length(t_s)]

source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)
nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))
strat1_w_GERD_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

for (i in 1:length(w)){
  strat1_w_f[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OCf
  strat1_w_GERD_f[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OCf
}


source('/Users/curtiu01/Dropbox/BaseCaseII/af_parameter_list.R')

femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)


parms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho =allfemales$rho,t=0)
parms['t'] = b
phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,b)
ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
surv80 <- ode$phi9[2]

surv=surv80

ts_star=69
parms['t'] = ts_star
phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,ts_star)
ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
surv_ts_star <- ode$phi9[2]
surv_ts_star_allfemales = surv_ts_star


BEonsetcum_MSCE80<-rep(1,b)
for (i in 1:(ts_star)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    parms['t'] = ts_star-xg[a]
    times = c(0,(ts_star-xg[a]))
    ode_i = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    GERDpr1=1-exp(-gerd1*min(gerd3,xg[a])-gerd2*max(0,xg[a]-gerd3))
    nu <- nu0*((1-GERDpr1)+RR*GERDpr1)
    GLquad1<-legauss(0 ,xg[a],20)
    wg2 <- GLquad1$weights
    xg2 <- GLquad1$mesh
    int1 <- rep(0,20)
    for ( j in 1:20){
      GERDpr=1-exp(-gerd1*min(gerd3,xg2[j])-gerd2*max(0,xg2[j]-gerd3))
      nu_temp <-nu0*((1-GERDpr)+RR*GERDpr)
      int1[j] <- nu_temp
    }
    nu_cum<-1-exp(-sum(wg2*int1))
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode_i$phi7[2]-ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcum_MSCE80[i] <- (1/(surv_ts_star-surv80))*(f_int)
}


## Compute SD probability 
SD_allfemales = BEonsetcum_MSCE80[ts_star]

### FOR FEMALES WITH GERD SYMPTOMS

parms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho =allfemales$rho,t=0)
parms['t'] = b
phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,b)
ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
surv80GERD <- ode$phi9[2]

ts_star=64
parms['t'] = ts_star
phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
times = c(0,ts_star)
ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
surv_ts_star <- ode$phi9[2]
surv_ts_starGERDfemales <-  surv_ts_star


BEonsetcumGERD_MSCE80<-rep(1,b)
for (i in 1:(ts_star)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    parms['t'] = ts_star-xg[a]
    times = c(0,(ts_star-xg[a]))
    ode_i = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    nu <- nu0*(RR*1)
    nu_cum<-1-exp(-nu*xg[a])
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode_i$phi7[2]-ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcumGERD_MSCE80[i] <- (1/(surv_ts_star-surv80GERD))*(f_int)
}

## Compute SD probability 
SD_GERDfemales=BEonsetcumGERD_MSCE80[ts_star]


################################################################
###    COMPUTE Overdiagnosis Function (OD(t_s^*)):           ###
##            Pr [ T_BE <= t_s | T_EAC > B]                  ###
################################################################

parms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho =allfemales$rho,t=0)

BEonsetcum_MSCE_noEAC80<-rep(0,b)
for (i in 1:(final_int)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  f4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    GERDpr1=1-exp(-gerd1*min(gerd3,xg[a])-gerd2*max(0,xg[a]-gerd3))
    nu <- nu0*((1-GERDpr1)+RR*GERDpr1)
    GLquad1<-legauss(0 ,xg[a],20)
    wg2 <- GLquad1$weights
    xg2 <- GLquad1$mesh
    int1 <- rep(0,20)
    for ( j in 1:20){
      GERDpr=1-exp(-gerd1*min(gerd3,xg2[j])-gerd2*max(0,xg2[j]-gerd3))
      nu_temp <-nu0*((1-GERDpr)+RR*GERDpr)
      int1[j] <- nu_temp
    }
    nu_cum<-1-exp(-sum(wg2*int1))
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4) 
  BEonsetcum_MSCE_noEAC80[i] <- (1/(surv80))*(f_int)
}

########################################
## Compute Over-screened probability  ##
########################################
ts_star = 69
prob_overdiag_noEAC80 = survdeath_by_OCf[(ts_star-9)]*(BEonsetcum_MSCE_noEAC80[ts_star]/strat1_w_f[(ts_star-9),5])*surv80
OD_allfemales = BEonsetcum_MSCE_noEAC80[ts_star]
os = prob_overdiag_noEAC80 

### GERD females only

BEonsetcumGERD_MSCE_noEAC80<-rep(0,b)
for (i in 1:(final_int)){
  GLquad<-legauss(0,i,gauss_points)
  xg <- GLquad$mesh
  wg <- GLquad$weights
  s4<-rep(0,gauss_points)
  for (a in 1:gauss_points){
    parms['t'] = b-xg[a]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,(b-xg[a]))
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    nu <- nu0*(RR*1)
    nu_cum<-1-exp(-nu*xg[a])
    nu_dens <- nu*(1-nu_cum)
    s4[a] <- (ode$phi7[2])*nu_dens
  }
  f_int<-sum(wg*s4)
  BEonsetcumGERD_MSCE_noEAC80[i] <- (1/(surv80GERD))*(f_int)
}
ts_star = 64
## Compute Over-screened probability 
prob_overdiagGERD_noEAC80 = survdeath_by_OCf[(ts_star-9)]*(BEonsetcumGERD_MSCE_noEAC80[ts_star]/strat1_w_GERD_f[(ts_star -9),5])*surv80GERD
OD_GERDfemales = BEonsetcumGERD_MSCE_noEAC80[ts_star]
osGERD = prob_overdiagGERD_noEAC80

## Compute PPV:


ss = 1 - os
ssGERD=1-osGERD
ppv_GERDfemale = ssGERD
ppv_female = ss

print(paste("PPV all females:",ppv_female))
print(paste("PPV GERD+ females:",ppv_GERDfemale))
print(paste("SD all females:",SD_allfemales))
print(paste("SD GERD+ females:",SD_GERDfemales))
print(paste("OD all females:",OD_allfemales))
print(paste("OD GERD+ females:",OD_GERDfemales))
print(paste("Cancer risk all females:",1-surv_ts_star_allfemales))
print(paste("Cancer risk GERD+ females:",1-surv_ts_starGERDfemales))


######################################################################
###       RESULTS FOR SENSITIVITY ANALYSES BORN 1950               ###
######################################################################
### Back to Strategy 1 w values
RR=5
b_c1 = 1950
start=10
w=seq(0,10, .25)
t_s=start:80
birth_cohort = b_c1
survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv[1:length(t_s)]
survdeath_by_OCf <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv[1:length(t_s)]


strat1_w_GERD_f_m = matrix(rep(0,1000*length(t_s)),ncol=length(t_s))
strat1_w_f_m = matrix(rep(0,1000*length(t_s)),ncol=length(t_s))
strat1_w_GERD_m = matrix(rep(0,1000*length(t_s)),ncol=length(t_s))
strat1_w_m = matrix(rep(0,1000*length(t_s)),ncol=length(t_s))

for (i in 1:1000){
  source('af_Bootstrap_sensitivity_2020.R')
  femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

  nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
  nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]


  p_EAC = (1-s_EAC(femaleparms,t_s))
  p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

  strat1_w_f_m[i,] = (nu_cum-w[5]*p_EAC)*survdeath_by_OCf
  strat1_w_GERD_f_m[i,] = (nu_cum_GERD-w[5]*p_EAC_GERD)*survdeath_by_OCf
  print(i)
  source('am_Bootstrap_sensitivity_2020.R')
  maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

  nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
  nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]


  p_EAC = (1-s_EAC(maleparms,t_s))
  p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

  strat1_w_m[i,] = (nu_cum-w[5]*p_EAC)*survdeath_by_OC
  strat1_w_GERD_m[i,] = (nu_cum_GERD-w[5]*p_EAC_GERD)*survdeath_by_OC
  print(i)
}
######################################################################
#### FOR TABLE 1, 95% confidence intervals on single age results #####
########  also used for Supplementary Figure S1A                 #####
######################################################################
quantile(t_s[apply(strat1_w_GERD_m,1,which.max)[1:1000]], prob=c(.025,.5,.975))
quantile(t_s[apply(strat1_w_m,1,which.max)[1:1000]], prob=c(.025,.5,.975))
quantile(t_s[apply(strat1_w_GERD_f_m,1,which.max)[1:1000]], prob=c(.025,.5,.975))
quantile(t_s[apply(strat1_w_f_m,1,which.max)[1:1000]], prob=c(.025,.5,.975))



##########################################################################################################################
##            FIGURE 3 for differing RR values (RR=2-6, baseline RR=5) in indepdendent validation                       ##
## versus published CORI data (Rubenstein et al 2010) + corresponding optimal ages from screening objective function.   ##
##########################################################################################################################

b_c1=1950
w=seq(0,10, .25)
t_s=start:80
birth_cohort = b_c1
survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv[1:length(t_s)]
survdeath_by_OCf <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv[1:length(t_s)]
RR=5


source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
  strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
  strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50=strat1_w_GERD_noMortality
strat1_w_noMortality50=strat1_w_noMortality

strat1_w_f_noMortality50=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50 = strat1_w_GERD_f_noMortality



RR=6
birth_cohort = 1950
source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}

for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
    strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50_6=strat1_w_GERD_noMortality
strat1_w_noMortality50_6=strat1_w_noMortality

strat1_w_f_noMortality50_6=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50_6 = strat1_w_GERD_f_noMortality



RR=2
birth_cohort = 1950
source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}

for (i in 1:length(w)){
    strat1_w_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OC
    strat1_w_GERD_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f_noMortality = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_f_noMortality[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts#*survdeath_by_OCf
    strat1_w_GERD_f_noMortality[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD#*survdeath_by_OCf
}


strat1_w_GERD_noMortality50_2=strat1_w_GERD_noMortality
strat1_w_noMortality50_2=strat1_w_noMortality

strat1_w_f_noMortality50_2=strat1_w_f_noMortality
strat1_w_GERD_f_noMortality50_2 = strat1_w_GERD_f_noMortality


## FIGURE 3 A-B
dev.new()
par(mfrow=c(1,2),xaxs="i", yaxs="i", las = 1, mar=c(5,5,2,3), bty='l')

plot(t_s[seq(16,66,1)],strat1_w_f_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="orange",lty=1,pch=18,cex=1,cex.axis=2, cex.lab=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="green3",lty=1,pch=17,cex=1)
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_f_noMortality50_2[seq(16,66),5]*100,strat1_w_f_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "orange", alpha.f = 0.1),border=adjustcolor( "orange", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_noMortality50_2[seq(16,66),5]*100,strat1_w_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "green3", alpha.f = 0.1),border=adjustcolor( "green3", alpha.f = 0.1))

legend('topleft',c('Model: Males',  'Model: Females' ), col=c("green3", "orange"),cex=1.5, lwd = 3, bty='n',lty=c(1))


plot(t_s[seq(16,66,1)],strat1_w_GERD_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="blue",lty=1,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "BE prevalence in non-EAC (%)", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_GERD_f_noMortality50[seq(16,66,1),5]*100,type='l',lwd=3,col="red",lty=1,pch=16,cex=1)
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_noMortality50_2[seq(16,66),5]*100,strat1_w_GERD_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "blue", alpha.f = 0.1),border=adjustcolor( "blue", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_f_noMortality50_2[seq(16,66),5]*100,strat1_w_GERD_f_noMortality50_6[seq(66,16),5]*100),col= adjustcolor( "red", alpha.f = 0.1),border=adjustcolor( "red", alpha.f = 0.1))

legend('topleft',c('Model: GERD Males',  'Model: GERD Females'), col=c("blue",  "red"),cex=1, lwd = 3, bty='n',lty=c(1))

## COMPARE WITH BEST3 for 69 year olds:

BE_prev_BEST3_compare =100*(strat1_w_GERD_noMortality50[seq(16,66,1),5][45]*.48+strat1_w_GERD_f_noMortality50[seq(16,66,1),5][45]*.52)
# [1] 6.285998

RR=5
source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}



for (i in 1:length(w)){
    strat1_w[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts*survdeath_by_OC
    strat1_w_GERD[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
  parms['t'] = t_s[j]
  phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
  times = c(0,t_s[j])
  ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
  surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
  strat1_w_f[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts*survdeath_by_OCf
  strat1_w_GERD_f[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD*survdeath_by_OCf
}


strat1_w_GERD_50=strat1_w_GERD
strat1_w_50=strat1_w

strat1_w_f_50=strat1_w_f
strat1_w_GERD_f_50 = strat1_w_GERD_f



RR=6
birth_cohort = 1950
source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts*survdeath_by_OC
    strat1_w_GERD[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}

for (i in 1:length(w)){
    strat1_w_f[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts*survdeath_by_OCf
    strat1_w_GERD_f[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD*survdeath_by_OCf
}


strat1_w_GERD_50_6=strat1_w_GERD
strat1_w_50_6=strat1_w

strat1_w_f_50_6=strat1_w_f
strat1_w_GERD_f_50_6 = strat1_w_GERD_f



RR=2
birth_cohort = 1950
source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
parms = maleparms

p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))

strat1_w_GERD = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allmales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}

for (i in 1:length(w)){
    strat1_w[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts*survdeath_by_OC
    strat1_w_GERD[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
parms= femaleparms

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))

strat1_w_GERD_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

surv_C_ts = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts[j]<- ode$phi9[2]
}
surv_C_ts_GERD = rep(0, length(t_s))
for (j in 1:length(t_s)){
    parms['t'] = t_s[j]
    phi0 = c(phi1=1, phi2=-allfemales$rho, phi3=1, phi4=0, phi5=1 , phi6=0, phi7=1,phi8=0, phi9=1, phi10=0)
    times = c(0,t_s[j])
    ode = as.data.frame(lsoda(phi0,times,func=phidot_GERD,parms,atol=1.e-9,rtol=1.e-9))
    surv_C_ts_GERD[j]<- ode$phi9[2]
}


for (i in 1:length(w)){
    strat1_w_f[,i] = (nu_cum-w[i]*p_EAC)/surv_C_ts*survdeath_by_OCf
    strat1_w_GERD_f[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)/surv_C_ts_GERD*survdeath_by_OCf
}


strat1_w_GERD_50_2=strat1_w_GERD
strat1_w_50_2=strat1_w

strat1_w_f_50_2=strat1_w_f
strat1_w_GERD_f_50_2 = strat1_w_GERD_f


## FIGURE 3 C-D
dev.new()
par(mfrow=c(1,2),xaxs="i", yaxs="i", las = 1, mar=c(5,5,2,3), bty='l')

plot(t_s[seq(16,66,1)],strat1_w_f_50[seq(16,66,1),5]*100,type='l',lwd=3,col="orange",lty=1,pch=18,cex=1,cex.axis=2, cex.lab=2, ylim=c(0,13),xlim=c(20,80),ylab= "Screening objective function", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_50[seq(16,66,1),5]*100,type='l',lwd=3,col="green3",lty=1,pch=17,cex=1)

polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_f_50_2[seq(16,66),5]*100,strat1_w_f_50_6[seq(66,16),5]*100),col= adjustcolor( "orange", alpha.f = 0.1),border=adjustcolor( "orange", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_50_2[seq(16,66),5]*100,strat1_w_50_6[seq(66,16),5]*100),col= adjustcolor( "green3", alpha.f = 0.1),border=adjustcolor( "green3", alpha.f = 0.1))

points(t_s[apply(strat1_w_50,2,which.max)[5]],strat1_w_50[apply(strat1_w_50,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)
points(t_s[apply(strat1_w_f_50,2,which.max)[5]],strat1_w_f_50[apply(strat1_w_f_50,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)

points(t_s[apply(strat1_w_50,2,which.max)[5]],0,lty=1,pch=18,cex=2, col="green3")
points(t_s[apply(strat1_w_f_50,2,which.max)[5]],0,col="orange",lty=1,pch=18,cex=2)

text(63.5,3,(bquote("t*"[s]~.("= 64",sep=''))), cex=2)
text(73,1.1,(bquote("t*"[s]~.("= 69",sep=''))), cex=2)

legend('topleft',c('Model: Males',  'Model: Females' ), col=c("green3", "orange"),cex=1.5, lwd = 3, bty='n',lty=c(1))

## GERD subpopulation
plot(t_s[seq(16,66,1)],strat1_w_GERD_50[seq(16,66,1),5]*100,type='l',lwd=3,col="blue",lty=1,pch=15,cex.axis=2, cex.lab=2, cex=2, ylim=c(0,13),xlim=c(20,80),ylab= "Screening objective function", xlab="Age")
lines(t_s[seq(16,66,1)],strat1_w_GERD_f_50[seq(16,66,1),5]*100,type='l',lwd=3,col="red",lty=1,pch=16,cex=1)

polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_50_2[seq(16,66),5]*100,strat1_w_GERD_50_6[seq(66,16),5]*100),col= adjustcolor( "blue", alpha.f = 0.1),border=adjustcolor( "blue", alpha.f = 0.1))
polygon(c(t_s[seq(16,66,1)],t_s[seq(66,16)]),c(strat1_w_GERD_f_50_2[seq(16,66),5]*100,strat1_w_GERD_f_50_6[seq(66,16),5]*100),col= adjustcolor( "red", alpha.f = 0.1),border=adjustcolor( "red", alpha.f = 0.1))
points(t_s[apply(strat1_w_GERD_50,2,which.max)[5]],strat1_w_GERD_50[apply(strat1_w_GERD_50,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)
points(t_s[apply(strat1_w_GERD_f_50,2,which.max)[5]],strat1_w_GERD_f_50[apply(strat1_w_GERD_f_50,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)
points(t_s[apply(strat1_w_GERD_50_2,2,which.max)[5]],strat1_w_GERD_50_2[apply(strat1_w_GERD_50_2,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)
points(t_s[apply(strat1_w_GERD_f_50_2,2,which.max)[5]],strat1_w_GERD_f_50_2[apply(strat1_w_GERD_f_50_2,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)
points(t_s[apply(strat1_w_GERD_50_6,2,which.max)[5]],strat1_w_GERD_50_6[apply(strat1_w_GERD_50_6,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)
points(t_s[apply(strat1_w_GERD_f_50_6,2,which.max)[5]],strat1_w_GERD_f_50_6[apply(strat1_w_GERD_f_50_6,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=3)

points(t_s[apply(strat1_w_GERD_50,2,which.max)[5]],0,col="blue",lty=1,pch=18,cex=2)
points(t_s[apply(strat1_w_GERD_f_50,2,which.max)[5]],0,col="red",lty=1,pch=18,cex=2, bg="green3")
text(56,7,(bquote("t*"[s]~.("= 58",sep=''))), cex=2)
text(63.5,3,(bquote("t*"[s]~.("= 64",sep=''))), cex=2)
legend('topleft',c('Model: GERD Males',  'Model: GERD Females'), col=c("blue",  "red"),cex=1, lwd = 3, bty='n',lty=c(1))





######################################
###     RE-SCREENING STRATEGY.     ###
######################################


############################################################################
###  RESULTS FOR MALES BORN 1950 - we only show these Results for RR=5.  ###
############################################################################

b_c1 = 1950
RR=5
start=10
birth_cohort = b_c1
ts_1 = c(45,50,58,64,69)
survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv

source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

i=1
ng=20
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp1 = survd*nu_cum2*(nu_cum1-comp2)
ts_21=ts_2

i=2
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp2 = survd*nu_cum2*(nu_cum1-comp2)
ts_22=ts_2

i=3
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp3 = survd*nu_cum2*(nu_cum1-comp2)
ts_23=ts_2

i=4
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp4 = survd*nu_cum2*(nu_cum1-comp2)
ts_24=ts_2

i=5
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp5 = survd*nu_cum2*(nu_cum1-comp2)
ts_25=ts_2

### GERD MALES ONLY
i=1
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp1_GERD = survd*nu_cum2*(nu_cum1-comp2)
ts_21_GERD=ts_2

i=2
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp2_GERD = survd*nu_cum2*(nu_cum1-comp2)
ts_22_GERD=ts_2

i=3
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp3GERD = survd*nu_cum2*(nu_cum1-comp2)
ts_23GERD=ts_2

i=4
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp4GERD = survd*nu_cum2*(nu_cum1-comp2)
ts_24GERD=ts_2

i=5
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OC[id2]-1+survdeath_by_OC[id1])/(survdeath_by_OC[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(maleparms,ts_2[j]-xg_tau))))
}

temp5GERD = survd*nu_cum2*(nu_cum1-comp2)
ts_25GERD=ts_2



###########################################################
### RESULTS FOR FEMALES BORN 1950 ###
###########################################################

b_c1 = 1950
ts_1 = c(45,50,58,64,69)

birth_cohort = b_c1
source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)
survdeath_by_OCf <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv


i=1
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp1f = survd*nu_cum2*(nu_cum1-comp2)
ts_21f=ts_2

i=2
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp2f = survd*nu_cum2*(nu_cum1-comp2)
ts_22f=ts_2

i=3
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp3f = survd*nu_cum2*(nu_cum1-comp2)
ts_23f=ts_2

i=4
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp4f = survd*nu_cum2*(nu_cum1-comp2)
ts_24f=ts_2

i=5
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp5f = survd*nu_cum2*(nu_cum1-comp2)
ts_25f=ts_2


### GERD FEMALES ONLY
i=1
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp1_GERDf = survd*nu_cum2*(nu_cum1-comp2)
ts_21_GERDf=ts_2

i=2
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp2_GERDf = survd*nu_cum2*(nu_cum1-comp2)
ts_22_GERDf=ts_2

i=3
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp3GERDf = survd*nu_cum2*(nu_cum1-comp2)
ts_23GERDf=ts_2

i=4
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp4GERDf = survd*nu_cum2*(nu_cum1-comp2)
ts_24GERDf=ts_2

i=5
ts_2=ts_1[i]:119
id1 = match(ts_1[i],start:119)
id2 = match(ts_2,start:119)
survd = (survdeath_by_OCf[id2]-1+survdeath_by_OCf[id1])/(survdeath_by_OCf[id1])
nu_cum1 = BEdensity_exactage_GERD(ts_2,gerd1,gerd2,gerd3,nu0,RR)$nu_cum - BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum 
nu_cum2=(1-BEdensity_exactage_GERD(ts_1[i],gerd1,gerd2,gerd3,nu0,RR)$nu_cum)
nu_cum = nu_cum1/nu_cum2
comp2 = rep(0,length(ts_2))
for (j in 1:length(comp2)){
  GLquadtau<-legauss(ts_1[i],ts_2[j],ng)
  wg_tau <- GLquadtau$weights
  xg_tau <- GLquadtau$mesh
  nu_dens = BEdensity_exactage_GERD(xg_tau,gerd1,gerd2,gerd3,nu0,RR)$nu_dens
  comp2[j]= sum(wg_tau*nu_dens*((1-s_MSCE(femaleparms,ts_2[j]-xg_tau))))
}

temp5GERDf = survd*nu_cum2*(nu_cum1-comp2)
ts_25GERDf=ts_2



print(paste("ts2 for ts1", ts_1[1],"all males:", ts_21[which.max(temp1)], "prob success:",max(temp1)))
print(paste("ts2 for ts1", ts_1[2],"all males:", ts_22[which.max(temp2)], "prob success:",max(temp2)))
print(paste("ts2 for ts1", ts_1[3],"all males:", ts_23[which.max(temp3)], "prob success:",max(temp3)))
print(paste("ts2 for ts1", ts_1[4],"all males:", ts_24[which.max(temp4)], "prob success:",max(temp4)))
print(paste("ts2 for ts1", ts_1[5],"all males:", ts_25[which.max(temp5)], "prob success:",max(temp5)))


print(paste("ts2 for ts1", ts_1[1],"GERD males:", ts_21_GERD[which.max(temp1_GERD)], "prob success:",max(temp1_GERD)))
print(paste("ts2 for ts1", ts_1[2],"GERD males:", ts_22_GERD[which.max(temp2_GERD)], "prob success:",max(temp2_GERD)))
print(paste("ts2 for ts1", ts_1[3],"GERD males:", ts_23GERD[which.max(temp3GERD)], "prob success:",max(temp3GERD)))
print(paste("ts2 for ts1", ts_1[4],"GERD males:", ts_24GERD[which.max(temp4GERD)], "prob success:",max(temp4GERD)))
print(paste("ts2 for ts1", ts_1[5],"GERD males:", ts_25GERD[which.max(temp5GERD)], "prob success:",max(temp5GERD)))

print(paste("ts2 for ts1", ts_1[1],"all females:", ts_21f[which.max(temp1f)], "prob success:",max(temp1f)))
print(paste("ts2 for ts1", ts_1[2],"all females:", ts_22f[which.max(temp2f)], "prob success:",max(temp2f)))
print(paste("ts2 for ts1", ts_1[3],"all females:", ts_23f[which.max(temp3f)], "prob success:",max(temp3f)))
print(paste("ts2 for ts1", ts_1[4],"all females:", ts_24f[which.max(temp4f)], "prob success:",max(temp4f)))
print(paste("ts2 for ts1", ts_1[5],"all females:", ts_25f[which.max(temp5f)], "prob success:",max(temp5f)))


print(paste("ts2 for ts1", ts_1[1],"GERD females:", ts_21_GERDf[which.max(temp1_GERDf)], "prob success:",max(temp1_GERDf)))
print(paste("ts2 for ts1", ts_1[2],"GERD females:", ts_22_GERDf[which.max(temp2_GERDf)], "prob success:",max(temp2_GERDf)))
print(paste("ts2 for ts1", ts_1[3],"GERD females:", ts_23GERDf[which.max(temp3GERDf)], "prob success:",max(temp3GERDf)))
print(paste("ts2 for ts1", ts_1[4],"GERD females:", ts_24GERDf[which.max(temp4GERDf)], "prob success:",max(temp4GERDf)))
print(paste("ts2 for ts1", ts_1[5],"GERD females:", ts_25GERDf[which.max(temp5GERDf)], "prob success:",max(temp5GERDf)))

###################################
###     Figure 4 (A-D)          ###
###################################

### GREYSCALE VERSION
dev.new() 
colors_o1 = gray.colors(5, start=0,end=0.7)
par(mfrow=c(2,2),mar=c(5,5,5,5))
#MALES
plot(0,xlim=c(0,119), ylim=c(0,(max(temp1)+.001)),col="white",xlab=expression("Screening time"~'t'['s2']~"(age)"),ylab="Re-screening objective function",yaxs="i",cex.lab=1.5,cex.axis=1.5)
lines(ts_21,temp1,col=colors_o1[1],lwd=2,lty=1)
lines(ts_22,temp2,col=colors_o1[2],lwd=2,lty=5)
lines(ts_23,temp3,col=colors_o1[3],lwd=2,lty=1)
lines(ts_24,temp4,col=colors_o1[4],lwd=2,lty=3)
lines(ts_25,temp5,col=colors_o1[5],lwd=2,lty=1)
legend("topleft",c(expression('t'['s1']==45),expression('t'['s1']==50),expression('t'['s1']==58),expression('t'['s1']==64),expression('t'['s1']==69)), col=colors_o1,lty=c(1,5,1,3,1),lwd=2,bty='n',cex=1)
points(ts_22[which.max(temp2)],max(temp2),col="black",cex=2,pch=18)
points(ts_21[which.max(temp1)],max(temp1),col="black",cex=2,pch=18)
points(ts_23[which.max(temp3)],max(temp3),col=1,cex=2,pch=18)
points(ts_24[which.max(temp4)],max(temp4),col=1,cex=2,pch=18)
points(ts_25[which.max(temp5)],max(temp5),col=1,cex=2,pch=18)

title(main="A)     All Males")

plot(0,xlim=c(0,119), ylim=c(0,(max(temp1_GERD)+.001)),col="white",xlab=expression("Screening time"~'t'['s2']~"(age)"),ylab="Re-screening objective function",yaxs="i",cex.lab=1.5,cex.axis=1.5)
lines(ts_21_GERD,temp1_GERD,col=colors_o1[1],lwd=2,lty=1)
lines(ts_22_GERD,temp2_GERD,col=colors_o1[2],lwd=2,lty=5)
lines(ts_23GERD,temp3GERD,col=colors_o1[3],lwd=2,lty=3)
lines(ts_24GERD,temp4GERD,col=colors_o1[4],lwd=2,lty=1)
lines(ts_25GERD,temp5GERD,col=colors_o1[5],lwd=2,lty=1)
legend("topleft",c(expression('t'['s1']==45),expression('t'['s1']==50),expression('t'['s1']==58),expression('t'['s1']==64),expression('t'['s1']==69)), col=colors_o1,lty=c(1,5,3,1,1),lwd=2,bty='n',cex=1.5)
points(ts_22_GERD[which.max(temp2_GERD)],max(temp2_GERD),col=1,cex=2,pch=18)
points(ts_21_GERD[which.max(temp1_GERD)],max(temp1_GERD),col=1,cex=2,pch=18)
points(ts_23GERD[which.max(temp3GERD)],max(temp3GERD),col=1,cex=2,pch=18)
points(ts_24GERD[which.max(temp4GERD)],max(temp4GERD),col=1,cex=2,pch=18)
points(ts_25GERD[which.max(temp5GERD)],max(temp5GERD),col=1,cex=2,pch=18)
title(main="B)     Males with GERD symptoms")

#FEMALES

plot(0,xlim=c(0,119), ylim=c(0,(max(temp1f)+.001)),col="white",xlab=expression("Screening time"~'t'['s2']~"(age)"),ylab="Re-screening objective function",yaxs="i",cex.lab=1.5,cex.axis=1.5)
lines(ts_21f,temp1f,col=colors_o1[1],lwd=2, lty=1)
lines(ts_22f,temp2f,col=colors_o1[2],lwd=2, lty=5)
lines(ts_23f,temp3f,col=colors_o1[3],lwd=2, lty =1)
lines(ts_24f,temp4f,col=colors_o1[4],lwd=2, lty=1)
lines(ts_25f,temp5f,col=colors_o1[5],lwd=2, lty=3)
legend("topleft",c(expression('t'['s1']==45),expression('t'['s1']==50),expression('t'['s1']==58),expression('t'['s1']==64),expression('t'['s1']==69)), col=colors_o1,lty=c(1,5,1,1,3),lwd=2,bty='n',cex=1.5)
points(ts_22f[which.max(temp2f)],max(temp2f),col=1,cex=2,pch=18)
points(ts_21f[which.max(temp1f)],max(temp1f),col=1,cex=2,pch=18)
points(ts_23f[which.max(temp3f)],max(temp3f),col=1,cex=2,pch=18)
points(ts_24f[which.max(temp4f)],max(temp4f),col=1,cex=2,pch=18)
points(ts_25f[which.max(temp5f)],max(temp5f),col=1,cex=2,pch=18)
title(main="C)     All Females")
plot(0,xlim=c(0,119), ylim=c(0,(max(temp1_GERDf)+.001)),col="white",xlab=expression("Screening time"~'t'['s2']~"(age)"),ylab="Re-screening objective function",yaxs="i",cex.lab=1.5,cex.axis=1.5)
lines(ts_21_GERDf,temp1_GERDf,col=colors_o1[1],lwd=2, lty=1)
lines(ts_22_GERDf,temp2_GERDf,col=colors_o1[2],lwd=2, lty=5)
lines(ts_23GERDf,temp3GERDf,col=colors_o1[3],lwd=2, lty=1)
lines(ts_24GERDf,temp4GERDf,col=colors_o1[4],lwd=2, lty=3)
lines(ts_25GERDf,temp5GERDf,col=colors_o1[5],lwd=2,lty=1)
legend("topleft",c(expression('t'['s1']==45),expression('t'['s1']==50),expression('t'['s1']==58),expression('t'['s1']==64),expression('t'['s1']==69)), col=colors_o1,lty=c(1,5,1,3,1),lwd=2,bty='n',cex=1.5)
points(ts_22_GERDf[which.max(temp2_GERDf)],max(temp2_GERDf),col=1,cex=2,pch=18)
points(ts_21_GERDf[which.max(temp1_GERDf)],max(temp1_GERDf),col=1,cex=2,pch=18)
points(ts_23GERDf[which.max(temp3GERDf)],max(temp3GERDf),col=1,cex=2,pch=18)
points(ts_24GERDf[which.max(temp4GERDf)],max(temp4GERDf),col=1,cex=2,pch=18)
points(ts_25GERDf[which.max(temp5GERDf)],max(temp5GERDf),col=1,cex=2,pch=18)
title(main="D)     Females with GERD symptoms")




######################################
###     SURVEILLANCE STRATEGY.     ###
######################################

#########################################################
### EAC RISK.  RESULTS FOR MALES BORN 1950            ###
#########################################################


tau = c(20, 30, 40)
ts_1 = c(45, 50, 58, 64)
t_F = 85
EACrisk=NULL
source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

## METHOD ONE ##
for (k in 1:length(tau)){
  EACrisk[[k]] = matrix(rep(0,length(ts_1)*length(1:t_F)),ncol=length(1:t_F))
  for ( i in 1:length(ts_1)){
    ts_2 = seq(ts_1[i],t_F,1)
    mean_p_stars = allmales$mu0*allmales$X*(ts_1[i]-tau[k])
    risk1 = 1-s_MSCE(maleparms,(ts_2-ts_1[i]))*s_3(maleparms,(ts_2-ts_1[i]))^mean_p_stars
    EACrisk[[k]][i,ts_1[i]:t_F]=risk1
  }
}

######################################
###  RESULTS FOR FEMALES BORN 1950 ###
######################################
ts_1 = c(45, 50, 64, 69)
EACriskf = NULL 
source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

## METHOD ONE ##
for (k in 1:length(tau)){
  EACriskf[[k]]=matrix(rep(0,length(ts_1)*length(1:t_F)),ncol=length(1:t_F))
  for ( i in 1:length(ts_1)){
    ts_2 = seq(ts_1[i],t_F,1)
    mean_p_stars = allfemales$mu0*allfemales$X*(ts_1[i]-tau[k])
    risk1 = 1-s_MSCE(femaleparms,(ts_2-ts_1[i]))*s_3(femaleparms,(ts_2-ts_1[i]))^mean_p_stars
    EACriskf[[k]][i,ts_1[i]:t_F]=risk1
  }
}

###################################
###     Figure 6 (A-B)          ###
###################################

dev.new()
par(mfrow=c(1,2), las=1)
colors_o1 = gray.colors(3, start=0,end=0.7)
ts_start=64
plot(ts_start:t_F,EACrisk[[1]][1,ts_start:t_F]*100,type="l",col="purple",lwd=1.5, xlab=expression("Screening time"~'t'['s2']~"(age)"),ylab="EAC Risk (%)",lty=1, cex.lab=1.25,cex.axis=1.5)
lines(ts_start:t_F,EACrisk[[2]][1,ts_start:t_F]*100,,col="forestgreen",lwd=1.5)
lines(ts_start:t_F,EACrisk[[3]][1,ts_start:t_F]*100,col="goldenrod",lwd=1.5)
lines(ts_start:t_F,EACrisk[[1]][2,ts_start:t_F]*100,type="l",col="purple",lwd=1.5,lty=2)
lines(ts_start:t_F,EACrisk[[2]][2,ts_start:t_F]*100,col="forestgreen",lwd=1.5,lty=2)
lines(ts_start:t_F,EACrisk[[3]][2,ts_start:t_F]*100,col="goldenrod",lwd=1.5,lty=2)
lines(ts_start:t_F,EACrisk[[1]][3,ts_start:t_F]*100,type="l",col="purple",lwd=1.5,lty=3)
lines(ts_start:t_F,EACrisk[[2]][3,ts_start:t_F]*100,col="forestgreen",lwd=1.5,lty=3)
lines(ts_start:t_F,EACrisk[[3]][3,ts_start:t_F]*100,col="goldenrod",lwd=1.5,lty=3)
lines(ts_start:t_F,EACrisk[[1]][4,ts_start:t_F]*100,type="l",col="purple",lwd=1.5,lty=4)
lines(ts_start:t_F,EACrisk[[2]][4,ts_start:t_F]*100,col="forestgreen",lwd=1.5,lty=4)
lines(ts_start:t_F,EACrisk[[3]][4,ts_start:t_F]*100,col="goldenrod",lwd=1.5,lty=4)
legend("topleft",c(expression('t'['s1']==45~','~tau ==20),expression('t'['s1']==45~','~tau ==30),expression('t'['s1']==45~','~tau ==40),expression('t'['s1']==50~','~tau ==20),expression('t'['s1']==50~','~tau ==30),expression('t'['s1']==50~','~tau ==40),expression('t'['s1']==58~','~tau ==20),expression('t'['s1']==58~','~tau ==30),
    expression('t'['s1']==58~','~tau ==40),expression('t'['s1']==64~','~tau ==20),expression('t'['s1']==64~','~tau ==30),expression('t'['s1']==64~','~tau ==40)), 
  col=c("purple","forestgreen","goldenrod","purple","forestgreen","goldenrod","purple","forestgreen","goldenrod"),lty=c(rep(1,3),2,2,2,3,3,3, rep(4,3)),lwd=1.5,bty='n',cex=1)
title(main="A)     All Males")
ts_start=69
b = t_F
plot(ts_start:b,EACriskf[[1]][1,ts_start:b]*100,type="l",col="purple",lwd=1.5, xlab=expression("Screening time"~'t'['s2']~"(age)"),ylab="EAC Risk (%)",lty=1,cex.lab=1.25,cex.axis=1.5)
lines(ts_start:b,EACriskf[[2]][1,ts_start:b]*100,,col="forestgreen",lwd=1.5)
lines(ts_start:b,EACriskf[[3]][1,ts_start:b]*100,col="goldenrod",lwd=1.5)
lines(ts_start:b,EACriskf[[1]][2,ts_start:b]*100,type="l",col="purple",lwd=1.5,lty=2)
lines(ts_start:b,EACriskf[[2]][2,ts_start:b]*100,col="forestgreen",lwd=1.5,lty=2)
lines(ts_start:b,EACriskf[[3]][2,ts_start:b]*100,col="goldenrod",lwd=1.5,lty=2)
lines(ts_start:b,EACriskf[[1]][3,ts_start:b]*100,type="l",col="purple",lwd=1.5,lty=3)
lines(ts_start:b,EACriskf[[2]][3,ts_start:b]*100,col="forestgreen",lwd=1.5,lty=3)
lines(ts_start:b,EACriskf[[3]][3,ts_start:b]*100,col="goldenrod",lwd=1.5,lty=3)
lines(ts_start:b,EACriskf[[1]][4,ts_start:b]*100,type="l",col="purple",lwd=1.5,lty=4)
lines(ts_start:b,EACriskf[[2]][4,ts_start:b]*100,col="forestgreen",lwd=1.5,lty=4)
lines(ts_start:b,EACriskf[[3]][4,ts_start:b]*100,col="goldenrod",lwd=1.5,lty=4)
legend("topleft",c(expression('t'['s1']==45~','~tau ==20),expression('t'['s1']==45~','~tau ==30),expression('t'['s1']==45~','~tau ==40),expression('t'['s1']==50~','~tau ==20),expression('t'['s1']==50~','~tau ==30),expression('t'['s1']==50~','~tau ==40),expression('t'['s1']==64~','~tau ==20),expression('t'['s1']==64~','~tau ==30),
    expression('t'['s1']==64~','~tau ==40),expression('t'['s1']==69~','~tau ==20),expression('t'['s1']==69~','~tau ==30),expression('t'['s1']==69~','~tau ==40)), 
  col=c("purple","forestgreen","goldenrod","purple","forestgreen","goldenrod","purple","forestgreen","goldenrod"),lty=c(rep(1,3),2,2,2,3,3,3, rep(4,3)),lwd=1.5,bty='n',cex=1)
title(main="B)     All Females")




##################################################################################################
###       RESULTS FOR Differences in predicted EAC risk based on t_s and tau above             ###
##################################################################################################
 
RR_EAC_f = EACriskf[[1]][2,85]/EACriskf[[3]][2,85]
RR_EAC_m = EACrisk[[1]][2,85]/EACrisk[[3]][2,85]

print(paste("Relative Risk (screen age 50) tau=20 to tau = 40, females:",RR_EAC_f))
print(paste("Relative Risk (screen age 50) tau=20 to tau = 40, males:",RR_EAC_m))



## Fitzerald et al BEST3 trial : 1750 patients underwent cytospong screening (1654 could swallow with a successful sample, 52% of whom where women), 
## Should condition to not include the 4 of whom had oesophago-gastric cancer, 127 total had BE
127/(1750-4)
## [1] 0.07273769
## Model estimate for same type of population
strat1_w_GERD_f_noMortality50[60,5]*0.52+strat1_w_GERD_noMortality50[60,5]*0.48
## [1] 0.06285998

## For Discussion
## Percent reduction in years of surveillance and in number of endoscopies if screening were to occur optimally
((0.67 - 0.61) / 0.67)*100
# [1] 8.955224
## Percent reduction in years of surveillance
((1.68 - 1.52) / 1.68)*100
# [1] 9.52381



######################################################################################################
###       RESULTS FOR SENSITIVITY ANALYSES BIRTH COHORTS 1960, 1970, 1980, 1990, 2000              ###
###                             Supplementary Figure s2                                            ###
######################################################################################################
RR=5
b_c1 = 1950
start=10
w=seq(0,10, .25)
t_s=start:80
birth_cohort = b_c1
survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv[1:length(t_s)]
survdeath_by_OCf <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv[1:length(t_s)]

source('/BEtoEAC_Results/am_parameter_list_2020.R')
maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)
nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]


p_EAC = (1-s_EAC(maleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))
strat1_w_GERD = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

for (i in 1:length(w)){
  strat1_w[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OC
  strat1_w_GERD[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OC
}


source('/BEtoEAC_Results/af_parameter_list_2020.R')
femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)
nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,nu0,RR)$nu_cum[t_s]
nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,nu0,RR)$nu_cum[t_s]

p_EAC = (1-s_EAC(femaleparms,t_s))
p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))
strat1_w_GERD_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
strat1_w_f = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

for (i in 1:length(w)){
  strat1_w_f[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OCf
  strat1_w_GERD_f[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OCf
}

birth_cohorts = seq(1960,2000,by=10)
RR=5
colors_bc = c('black','brown','grey','purple','cyan')
dev.new()

par(xaxs="i", yaxs="i", las = 1, mar=c(5,5,5,5))
plot(t_s[seq(15,65,1)],strat1_w_GERD[seq(15,65,1),5]*100,type='l',lwd=4,col="blue",lty=1,pch=15,cex.axis=2, cex.lab=2, cex=1, ylim=c(0,12),xlim=c(20,80),ylab= "Optimal BE yield (%)", xlab="Age", main= "Optimal Screen Ages")
lines(t_s[seq(15,65,1)],strat1_w_GERD_f[seq(15,65,1),5]*100,type='l',lwd=4,col="red",lty=1,pch=16,cex=1)


lines(t_s[seq(15,65,1)],strat1_w_f[seq(15,65,1),5]*100,type='l',lwd=4,col="orange",lty=1,pch=18,cex=1)
lines(t_s[seq(15,65,1)],strat1_w[seq(15,65,1),5]*100,type='l',lwd=4,col="green3",lty=1,pch=17,cex=1)


points(t_s[apply(strat1_w_GERD,2,which.max)[5]],strat1_w_GERD[apply(strat1_w_GERD,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=5)
points(t_s[apply(strat1_w_GERD_f,2,which.max)[5]],strat1_w_GERD_f[apply(strat1_w_GERD_f,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=5)

points(t_s[apply(strat1_w,2,which.max)[5]],strat1_w[apply(strat1_w,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=5)
points(t_s[apply(strat1_w_f,2,which.max)[5]],strat1_w_f[apply(strat1_w_f,2,which.max)[5],5]*100,col="green",lty=1,pch=18,cex=5)

points(t_s[apply(strat1_w_GERD,2,which.max)[5]],0,col="blue",lty=1,pch=18,cex=2)
points(t_s[apply(strat1_w,2,which.max)[5]],0,col="green3",lty=1,pch=18,cex=2)
points(t_s[apply(strat1_w_f,2,which.max)[5]],0,col="orange",lty=1,pch=18,cex=2)

legend('topleft',c('Male GERD','Female GERD','Males','Females'), col=c("blue", "red", "green3",  "orange"),pch=c(15,16,17,18),cex=1.5, bty='n',lty=c(1))

text(56,7,"t* = 58", cex=2)
text(63.5,3,"t* = 64", cex=2)
text(73,1.1,"t* = 69", cex=2)



for (k in 1:length(birth_cohorts)){

    b_c1 = birth_cohorts[k]
    birth_cohort = b_c1
    source('/BEtoEAC_Results/am_parameter_list_2020.R')
    maleparms = c(nu0=allmales$nu, alphap=allmales$alphaP,betap=allmales$betaP,alpham=allmales$alphaM,betam=allmales$betaM,mu0= X*allmales$mu0,mu1=allmales$mu1,mu2=allmales$mu2,rho=allmales$rho,t=0)

    nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0,RR)$nu_cum[t_s]
    nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,amnu0, RR)$nu_cum[t_s]


    p_EAC = (1-s_EAC(maleparms,t_s))
    p_EAC_GERD = (1-s_EAC_GERD(maleparms,t_s))
    #survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv
    survdeath_by_OC <- AllMortalityGenUS(race="all",byear=birth_cohort, start=start)$surv[1:length(t_s)]

    strat1_w_GERDtemp = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
    strat1_wtemp = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

    for (i in 1:length(w)){
      strat1_wtemp[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OC
      strat1_w_GERDtemp[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OC
    }

    start=10
    w=seq(0,20, .25)
    t_s=start:119
    birth_cohort = b_c1

    source('/BEtoEAC_Results/af_parameter_list_2020.R')
    femaleparms = c(nu0=allfemales$nu, alphap=allfemales$alphaP,betap=allfemales$betaP,alpham=allfemales$alphaM,betam=allfemales$betaM,mu0= X*allfemales$mu0,mu1=allfemales$mu1,mu2=allfemales$mu2,rho=allfemales$rho,t=0)

    nu_cum = BEdensity(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]
    nu_cum_GERD = BEdensity_GERD(t_s[length(t_s)],gerd1,gerd2,gerd3,afnu0,RR)$nu_cum[t_s]

    p_EAC = (1-s_EAC(femaleparms,t_s))
    p_EAC_GERD = (1-s_EAC_GERD(femaleparms,t_s))
    survdeath_by_OC <- AllMortalityGenUS(race="all",sex="female",byear=birth_cohort, start=start)$surv[1:length(t_s)]
    strat1_w_GERD_ftemp = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))
    strat1_w_ftemp = matrix(rep(0,length(w)*length(t_s)),ncol=length(w))

    for (i in 1:length(w)){
      strat1_w_ftemp[,i] = (nu_cum-w[i]*p_EAC)*survdeath_by_OC
      strat1_w_GERD_ftemp[,i] = (nu_cum_GERD-w[i]*p_EAC_GERD)*survdeath_by_OC
    }

    lines(t_s[seq(15,65,1)],strat1_w_GERDtemp[seq(15,65,1),5]*100,type='l',lwd=3,col=colors_bc[k],lty=(1+k),cex.axis=2, cex.lab=2, cex=2, ylim=c(0,12),xlim=c(20,80),ylab= "Optimal BE yield (%)", xlab="Age", main= "Optimal Screen Ages")
    lines(t_s[seq(15,65,1)],strat1_w_GERD_ftemp[seq(15,65,1),5]*100,type='l',lwd=3,col=colors_bc[k],lty=(1+k),cex=1)


    lines(t_s[seq(15,65,1)],strat1_w_ftemp[seq(15,65,1),5]*100,type='l',lwd=3,col=colors_bc[k],lty=(1+k),cex=1)
    lines(t_s[seq(15,65,1)],strat1_wtemp[seq(15,65,1),5]*100,type='l',lwd=3,col=colors_bc[k],lty=(1+k),cex=1)


    points(t_s[apply(strat1_w_GERDtemp,2,which.max)[5]],strat1_w_GERDtemp[apply(strat1_w_GERDtemp,2,which.max)[5],5]*100,col=colors_bc[k],lty=(1+k),pch=(14+k),cex=3)
    print(paste("GERD males:", t_s[apply(strat1_w_GERDtemp,2,which.max)[5]], "prob success: ",strat1_w_GERDtemp[apply(strat1_w_GERDtemp,2,which.max)[5],5]*100, "birth cohort: ", b_c1))
    
    points(t_s[apply(strat1_w_GERD_ftemp,2,which.max)[5]],strat1_w_GERD_ftemp[apply(strat1_w_GERD_ftemp,2,which.max)[5],5]*100,col=colors_bc[k],lty=(1+k),pch=(14+k),cex=3)
    print(paste("GERD females:", t_s[apply(strat1_w_GERD_ftemp,2,which.max)[5]], "prob success: ",strat1_w_GERD_ftemp[apply(strat1_w_GERD_ftemp,2,which.max)[5],5]*100, "birth cohort: ", b_c1))
    
    points(t_s[apply(strat1_wtemp,2,which.max)[5]],strat1_wtemp[apply(strat1_wtemp,2,which.max)[5],5]*100,col=colors_bc[k],lty=(1+k),pch=(14+k),cex=3)
    print(paste("males:", t_s[apply(strat1_wtemp,2,which.max)[5]], "prob success: ",strat1_wtemp[apply(strat1_wtemp,2,which.max)[5],5]*100, "birth cohort", b_c1))

    points(t_s[apply(strat1_w_ftemp,2,which.max)[5]],strat1_w_ftemp[apply(strat1_w_ftemp,2,which.max)[5],5]*100,col=colors_bc[k],lty=(1+k),pch=(14+k),cex=3)
    print(paste("females:", t_s[apply(strat1_w_ftemp,2,which.max)[5]], "prob success: ",strat1_w_ftemp[apply(strat1_w_ftemp,2,which.max)[5],5]*100, "birth cohort: ", b_c1))
}

legend('topright',as.character(birth_cohorts),col=colors_bc,pch=15:19,lty=2:6, lwd=3,bty='n',cex=1.5)

