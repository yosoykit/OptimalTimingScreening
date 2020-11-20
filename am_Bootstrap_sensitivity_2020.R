############################################################################################## 
### BOOTSRAP ANALYSIS OF PARAMETERS SAMPLED FROM MEAN AND 95% INTERVALS OF MCMC ESTIMATES) ###
###									 USE VALUES FROM 2014 CEBP paper                       ###
##############################################################################################


######################################################################################################
##  EXAMPLE PARAMETERS PROVIDED. Utilized for Screening paper 2014, fit to SEER EAC incidence data ###
##                			  all MALE PARAMETERS: AGE-COHORT MODEL									##
##  			birth cohort variables for sigmoidal g_P and g_M: See Supplementary 				##
######################################################################################################

param_list <- function(X, nu, alphaP, alphaM, rho, gP,gM,mu0,mu1,mu2){
	betaM <- alphaM-gM-rho
	mu2eff <- mu2*(1-betaM/alphaM)				## effective rate to pre-clinical cancer 
	betaeff<- alphaP-gP-mu2eff
	betaP<- alphaP-gP-mu2
	peff <- (1/2)*(-alphaP+betaeff + mu2eff - sqrt((alphaP+betaeff+mu2eff)^2-4*alphaP*betaeff))
	qeff <- (1/2)*(-alphaP+betaeff + mu2eff + sqrt((alphaP+betaeff+mu2eff)^2-4*alphaP*betaeff))
	pP <- (1/2)*(-alphaP+betaP + mu2 - sqrt((alphaP+betaP+mu2)^2-4*alphaP*betaP))
	qP <- (1/2)*(-alphaP+betaP + mu2 + sqrt((alphaP+betaP+mu2)^2-4*alphaP*betaP))
	pM <- (1/2)*(-alphaM+betaM + rho - sqrt((alphaM+betaM+rho)^2-4*alphaM*betaM))
	qM <- (1/2)*(-alphaM+betaM + rho + sqrt((alphaM+betaM+rho)^2-4*alphaM*betaM))
	zeta=-qM/pM 
	tlag <- log((zeta/(1+zeta)))/pM     				## exact T_2
	return(list(X=X,nu=nu,mu0=mu0,mu1=mu1,mu2=mu2, rho=rho,alphaP=alphaP,alphaM=alphaM, betaP=betaP,betaM=betaM, gP=gP,gM=gM, tlag=tlag, mu2eff=mu2eff,betaeff=betaeff,pP=pP, qP=qP,pM=pM,qM=qM,peff=peff,qeff=qeff))
}



 

##########################################################################################################
#### PARAMETERS BELOW WERE ESTIMATED WITH MLE (+ MCMC 95% CI) in TABLE 1 SCREENING PAPER 
##########################################################################################################
amnu0<- 0.36494E-03 
amnu0SD = (0.000413-0.000319)/3.92
amnu0=rnorm(1,amnu0, amnu0SD)
if (amnu0<=0){
	while(amnu0<=0){
		amnu0<- 0.36494E-03 
		amnu0=rnorm(1,amnu0, amnu0SD)
	}
}

ammu0<-  0.79942E-03
ammu0SD <- (0.000983-0.000638)/3.92
ammu0temp = rnorm(1,ammu0, ammu0SD)
if (ammu0temp<=0){
	while(ammu0temp<=0){
		ammu0temp = rnorm(1,ammu0, ammu0SD)
	}
}
crypts_per_avBE<- (1+15*(.5/4.5))*10*75/15*1000
ammu0<- ammu0temp*250000/crypts_per_avBE
ammu1 <-  ammu0temp


ammu2<-    0.45439E-04 
ammu2SD = (0.0000647-0.0000365)/3.92
ammu2=rnorm(1,ammu2, ammu2SD)
if (ammu2<=0){
	while(ammu2<=0){
		ammu2<-    0.45439E-04 
		ammu2=rnorm(1,ammu2, ammu2SD)
	}
}

am_g0<- 0.99061E-01      
am_g0SD <-  (1.099E-01 -0.928E-01)/3.92
am_g0 <- rnorm(1, am_g0, am_g0SD)
if (am_g0<=0){
	while(am_g0<=0){
		am_g0<- 0.99061E-01      
		am_g0 <- rnorm(1, am_g0, am_g0SD)
	}
}

am_g1<- 0.50880E+00     
am_g1SD <-  (0.59E+00 -0.275E+00)/3.92
am_g1 <- rnorm(1, am_g1, am_g1SD)
if (am_g1<=0){
	while(am_g1<=0){
		am_g1<- 0.50880E+00     
		am_g1 <- rnorm(1, am_g1, am_g1SD)
	}
}

am_g2<- 0.53768E-01  
am_g2SD <-  (0.572E-01 -0.483E-01)/3.92
am_g2 <- rnorm(1, am_g2, am_g2SD)
if (am_g2<=0){
	while(am_g2<=0){
	am_g2<- 0.53768E-01  
	am_g2 <- rnorm(1, am_g2, am_g2SD)
	}
}

am_trefc<- 0.19125E+04 
am_trefcSD <- (1914.1-1909.1)/3.92
am_trefc <- rnorm(1,am_trefc,am_trefcSD)
if (am_trefc<=0){
	while(am_trefc<=0){
		am_trefc<- 0.19125E+04 
		am_trefc <- rnorm(1,am_trefc,am_trefcSD)
	}
}


########################################################
am_gC0<- 0.75000E+00 
kstem<- 4    				## stem cells/crypt
#X<- 250000*kstem 			## total stem cells in a 5 cm BE segment assuming 250,000 crypts/5cm
X <- crypts_per_avBE*kstem

amrho<- 1.0000E-09       	## detection rate of clinical cancers
amalphaP0<- 10				## premalignant cell division rate
amalphaM0 <- 150			## pre-clinical malignant cell division rate
birth_cohort<- 1950
gerd1<- 0.00061422    
gerd2<- 0.0070447     
gerd3<- 26.002   


am_gP<- am_g0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
am_gM<- am_gC0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
#am_gM=am_gM*1.25
amalphaP<- amalphaP0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
amalphaM<- amalphaM0*(am_g1+2/(1+exp(-am_g2*(birth_cohort-am_trefc))))
allmales<- param_list(X, amnu0, amalphaP, amalphaM, amrho, am_gP,am_gM, ammu0, ammu1, ammu2)
amparams<- c(allmales$mu0,allmales$mu1,allmales$betaP,allmales$mu2,allmales$alphaP,allmales$betaM,allmales$rho,allmales$alphaM,allmales$pM,allmales$qM,allmales$tlag)

