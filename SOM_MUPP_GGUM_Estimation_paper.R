#Supplemental Online Material: 
#R Code for the MML-EM estimation of the MUPP-GGUM  
 
#set the working directory where the "SOM_MUPP_GGUM_Estimation_paper.Rdata." file is saved, for example, "c:\MUPP_GGUM"
setwd("c:/MUPP_GGUM") 

#load data (workspace) and library  
load("SOM_MUPP_GGUM_Estimation_paper.Rdata")

library(mirt)
library(mvtnorm)

#In the code below, the functions I wrote are explicitly specified or noted; all the other functions are from the mirt package or R. 
#You can always use the "dump" function to see the source code of a function. For example,

#dump("mirt", file="", control= NULL)

#############################
#simulate MUPP-GGUM data 
#############################

#number of traits
nf <- 3

#number of statements measuring each trait
nstate.per.f <- 20

#number of items
nitem <- nf*nstate.per.f/2

#sample size
nsample <-3000

#correlation matrix of trait scores
theta.cor <- matrix(0.3, nr=nf, nc=nf)
diag(theta.cor) <- 1

#A nitem X two statements matrix indicating the trait measured by each statement in each item
dim <- matrix(rep(combn(nf,2), nstate.per.f/2), nc=2, byrow=T)

##################################################################
#Function to create lognormal random variables by defining mean and variance at the normal scale 
#
#Description:
#Generate lognormal random variables with mean and variance input at the normal scale 
#
#Usage:
#lognorm(N, M, V)
#
#Arguments:
#N -- number of random variables to be drawn
#M -- mean of the lognormal distribution at the normal scale
#V -- variance of the lognormal distribution at the normal scale
#
#Return: a numeric vector with N random variables
##################################################################

#draw discrimination parameters from a lognormal distribution with mean=2 and sd=0.2 at the normal scale
a <- round(lognorm(nf*nstate.per.f, 2, 0.2^2),2)
a1 <- matrix(a, nc=2)

#draw location parameters from a uniform distribution between 0.5 and 1.5; change the second statements' location 
#parameters to negative
d <- round(runif(nf*nstate.per.f, 0.5, 1.5),2)
d1 <- matrix(d, nc=2)
d1[,2] <- -d1[,2]

#draw step parameters from a normal distribution with mean=1 and sd=1
#set the first statements' step parameters to 0
t <- round(rnorm(nf*nstate.per.f/2, 1, 1),2)
t1 <- cbind(rep(0,length(t)),t)

#combine all item parameters to a nitem X 6 matrix
#the columns in order are the first statement's discrimination, location, step, and
#the second statement's
par <- cbind(a1,d1, t1)[,c(1,3,5,2,4,6)]

#draw trait scores from a standard multivariate normal with 
#an intertrait correlation matrix = theta.cor
theta <- matrix(round(rmvnorm(nsample,rep(0, nf),theta.cor),2), ncol=nf)

##################################################################
#MUPP-GUMM's item response function (Equation 1 in the paper) 
#
#Description:
#Calculate an item's response probabilities based on MUPP-GGUM  
#
#Usage:
#P.MUPP.ggum.s(par,theta1, theta2)
#
#Arguments:
#par -- a vector including 6 item parameter in an item in order:
#	the first statement's discrimination, location, step, and
#	the second statement's
#theta1 -- a vector with nsample scores of the trait measured by the first statement
#theta2 -- a vector with nsample scores of the trait measured by the second statement
#
#Return: a nsample X 2 matrix containing the response probabilities for scores 0 and 1, respectively
##################################################################

#Generate data (a nsample X nitem matrix containing item responses) based on 
#MUPP-GGUM; a score category of an item must have at least 5 responses
simdata <- matrix(NA, nr=nsample, nc=nrow(par))
repeat{
for (i1 in 1:nrow(par)){
		prob <- P.MUPP.ggum.s(par[i1,], theta[,dim[i1,1]], theta[,dim[i1,2]])
		simdata[,i1] <- rbinom(n=nsample, size=1, prob[,2])		
}
if (!any(apply(simdata, 2, table)<5)) break
}
colnames(simdata) <- paste0("item", 1:nrow(par))

#Create new item type "MUPP_ggum" in mirt
name <- 'MUPP_ggum'
para <- rep(0, 3*nf)

#a=alpha, b=delta, t=tao
names(para)<- c(paste0(c("a","b","t"), rep(1:nf, each=3)))
est <- rep(F, 3*nf)

#P.MUPP.ggum is a wrapper of the MUPP-GGUM's item response function (Equation 1 in the paper).
#P.MUPP.ggum.gr is a wrapper of the first derivative function of the MUPP-GGUM's item response function, from the symbolic derivative function, Deriv.
#P.MUPP.ggum.hss is a wrapper of the second derivative function of the MUPP-GGUM's item response function, from the symbolic derivative function, Deriv.
Item.MUPP.ggum <- createItem(name, par=para, est=est, P=P.MUPP.ggum, gr=P.MUPP.ggum.gr, hss=P.MUPP.ggum.hss)

#create the model object in mirt  
Q <- matrix(0, nr=nitem, ncol=3, dimnames = list(NULL, paste0('F', 1:3)))
for (i in 1:nrow(Q)){
	Q[i,dim[i,]] <-1
}
Model.MUPP.sim <- mirt.model(Q)

#create user data used in mirt
user.data <- vector('list', nitem)
for ( i in 1:nitem){
	user.data[[i]] <- dim[i,]
}

#combine all item parameters and statements' trait IDs to a nitem X 8 matrix
#the columns in order are the first statement's discrimination, location, step, trait ID, and
#the second statement's
par <- cbind(par, dim)[,c(1,2,3,7,4,5,6,8)]
colnames(par) <- paste0(c("a", "b", "t","dim"), rep(1:2, each=4))

#output the default starting values of item parameters from mirt
sv <- mirt(simdata, Model.MUPP.sim, 'MUPP_ggum', customItems=list(MUPP_ggum=Item.MUPP.ggum), pars = 'values')

#set the starting values of item parameters to the true values and fix the intertrait correlations to the true values
tmp1 <- NULL
j <- 1
for ( i in 1:nrow(par)){		
	p <-""
	tmp1 <- rbind(tmp1, data.frame(item=paste0("item", i, p), name=paste0(c("a"), par[j,"dim1"]), value= unlist(par[j,c("a1")]),est=c(T),lbound=0))
	tmp1 <- rbind(tmp1, data.frame(item=paste0("item", i, p), name=paste0(c("a"), par[j,"dim2"]), value= unlist(par[j,c("a2")]),est=c(T),lbound=0))
	tmp1 <- rbind(tmp1, data.frame(item=paste0("item", i, p), name=paste0(c("b"), par[j,"dim1"]), value= unlist(par[j,c("b1")]),est=c(T),lbound=-Inf))
	tmp1 <- rbind(tmp1, data.frame(item=paste0("item", i, p), name=paste0(c("b"), par[j,"dim2"]), value= unlist(par[j,c("b2")]),est=c(T),lbound=-Inf))
	tmp1 <- rbind(tmp1, data.frame(item=paste0("item", i, p), name=paste0(c("t"), par[j,"dim1"]), value= unlist(par[j,c("t1")]),est=c(F),lbound=-Inf))
	tmp1 <- rbind(tmp1, data.frame(item=paste0("item", i, p), name=paste0(c("t"), par[j,"dim2"]), value= unlist(par[j,c("t2")]),est=c(T),lbound=-Inf))		
	j <- j+1
}
sv2 <- sv
for (i in 1:nrow(tmp1)){
	sv2[sv2[,"item"]==tmp1[i,"item"] & sv2[,"name"]==tmp1[i,"name"], c("value", "est", "lbound")] <- tmp1[i,c("value", "est", "lbound")]
}
tmp3 <- round(theta.cor,2)
sv2[grepl("COV", sv2$name, fixed=T), ]$value <- tmp3[lower.tri(tmp3,T)]
sv2.MUPP.ggum.sim <- sv2

#head(sv2,20)
#tail(sv2,50)

#Estimate MUPP-GGUM using MML-EM method
sim.ggum.MUPP.est <- mirt(simdata, Model.MUPP.sim, 'MUPP_ggum', customItems=list(MUPP_ggum=Item.MUPP.ggum), pars = sv2.MUPP.ggum.sim, 
  large=T, customItemsData=user.data, SE=F, TOL=0.001) 

#retrieve item parameter estimates 
tmp <- coef(sim.ggum.MUPP.est, simplify=T)[[1]]

tmp1 <- NULL
for (i in 1:nrow(tmp)){
	tmp1 <- rbind(tmp1, tmp[i, paste0(c("a","b", "t"), rep(dim[i,],each=3))])
}
par.est <- as.matrix(tmp1)

#retrieve estimation time 
extract.mirt(sim.ggum.MUPP.est, "time")

#retrieve convergence status 
extract.mirt(sim.ggum.MUPP.est,"converged")

#Estimate MAP trait scores 
score.est <- fscores(sim.ggum.MUPP.est, method="MAP", full.scores.SE=F)


