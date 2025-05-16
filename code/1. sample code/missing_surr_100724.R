
########################################
### Missing Surrogate Idea #############
########################################

# To get things started, let's just do 1 setting where 1) the true proportion explained (R) is 0.50, and 2) the data are generated from a parametric model. We can move on to other settings after this one is working. Let's call this setting 1. 
#this maps to this paper, DOI: 10.1002/sim.6820 setting 1 described in Sectiin 5.1, with parameter values that make R=0.50. It just takes some algebra to show the truth is 0.50.

########################################
### Functions to generate data ######
########################################

gen.data = function(setting, n1, n0) {
	if(setting == 1)
	{	s1 = g.1(n1)
		y1 = f.cond.1(s1)
		s0 = g.0(n0)
		y0 = f.cond.0(s0)
		return(list("s1" = s1, "y1" = y1, "s0" = s0, "y0" = y0))
	}
	if(setting == 2)
	{
	s1 = g.1.exp(n1)
	y1 = f.cond.1.exp(s1)
	s0 = g.0.exp(n0)
	y0 = f.cond.0.exp(s0)
	return(list("s1" = s1, "y1" = y1, "s0" = s0, "y0" = y0))
	}
	if(setting==3)
	{
	s1 = g.1(n1, alpha0=-0.33333333333)
	y1 = f.cond.1(s1)
	s0 = g.0(n0, alpha0=-0.33333333333)
	y0 = f.cond.0(s0)
	return(list("s1" = s1, "y1" = y1, "s0" = s0, "y0" = y0))
	}

}
	
f.cond.1 = function(s.vector) {
	eps1 = rnorm(length(s.vector),0,3)
	y1 = 2+5*s.vector+1 + 1*s.vector + eps1
	return(y1)		
}

f.cond.0 = function(s.vector) {
	eps0 = rnorm(length(s.vector),0,3)
	y0 = 2+5*s.vector+ eps0
	return(y0)		
}

f.cond.0.exp = function(s.vector ) {
	y = 0.8*s.vector^2 + .2*exp(s.vector/5) + exp(rnorm(length(s.vector),0 ,.3))
	return(y)
}

f.cond.1.exp = function(s.vector) {
	 y = .5+s.vector^2 + .5*exp(s.vector/5) + exp(rnorm(length(s.vector),0, .3))
	return(y)		
}

g.1 = function(n, alpha0=5) { return(rnorm(n, alpha0 + 1,2))}
g.0 = function(n,alpha0=5) { return(rnorm(n, alpha0,1))}
g.1.exp = function(n) { return(exp(rnorm(n,1.7, .2)))}
g.0.exp = function(n) { return(exp(rnorm(n,1.62, .1)))}

########################################
### Function to generate the truth ####
########################################
#this function will get more complicated when the parametric models don't hold

get.true = function(setting) {
	if(setting == 1)
	{return(list("delta.s.true" = 6, "delta.true" = 12, "R.true" = 0.5))}	

	if(setting == 2)
	{data = gen.data(setting=2, n1=100000, n0=100000) 
	s1 = data$s1
	y1 = data$y1
	s0 = data$s0
	y0 = data$y0
	delta.true = mean(y1)-mean(y0)
	y1.ifs0 = f.cond.1.exp(s0)
	delta.s.true = mean(y1.ifs0) - mean(y0)
	R.true = 1-(delta.s.true/delta.true)
	return(list("delta.s.true" = delta.s.true, "delta.true" = delta.true, "R.true" = R.true))
	}
	
	if(setting == 3)
	{return(list("delta.s.true" = .66666666666666666, "delta.true" = 6.666666666666666, "R.true" = 0.9))}	

}


########################################
### Kernel function ####
########################################
# it just takes an x and pops back the kernel

kf=function(x, h){return(dnorm(x/h)/h)}

########################################
### Some examples/explanation ####
########################################
#x is measuring how far away point a is from b, h is the bandwidth. Let's say the bandwidth is 0.2. And let's say 1 is 0.1 away from b, then the weight that a gets is:
kf(0.1,0.2)
#now what if a is 0.9 away from b, that is much farther, the weight it gets is: 
kf(0.9,0.2)
#much smaller because it's so far away


######################################################################
### Kernel function but setting it up to deal with matrix/vector ####
######################################################################
#but same idea as kf above

VTM<-function(vc, dm){
	#takes vc and makes it the repeated row of a matrix, repeats it dm times
     matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
    }
    
Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix ##
  { 
    out = (VTM(zz,length(zi))- zi)/bw
    return(dnorm(out)/bw)
           
  }

########################################
### My functions to estimate delta.s ####
########################################

#I make this so I can use an apply function in the next function
pred.smooth <-function(zz,zi.one, bw=NULL,y1, weight=NULL) {
  	if(is.null(bw)) { bw = bw.nrd(zz)/((length(zz))^(.10))}
  	if(is.null(weight)) {weight = rep(1, length(y1))}
  	return(sum(Kern.FUN(zz,zi.one,bw=bw)*y1*(1/weight))/sum(Kern.FUN(zz,zi.one,bw=bw)*1/weight))
}


delta.s.single = function(sone,szero,yone,yzero, h.select = NULL, weight.1 = NULL, weight.0=NULL, n0.all=NULL) {
	#we can talk about the bandwidth later, this default should work ok
	if(is.null(h.select)) {h.select = bw.nrd(sone)*(length(sone)^(-0.25)) 
	}
	if(is.null(weight.1) & is.null(weight.0)){
		mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone)
			if(sum(is.na(mu.1.s0))>0){
  				print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
    		c.mat = cbind(szero, mu.1.s0)
    		for(o in 1:length(mu.1.s0)) {
    			if(is.na(mu.1.s0[o])){
    				distance = abs(s0.new - s0.new[o])
    				c.temp = cbind(c.mat, distance)
    				c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    				new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    				mu.1.s0[o] = new.est[1]   #in case there are multiple matches
    				}
}}
		delta.s = mean(mu.1.s0) - mean(yzero)
	}
	if(!is.null(weight.1) & !is.null(weight.0)){
		mu.1.s0 = sapply(szero,pred.smooth,zz=sone, bw=h.select, y1=yone, weight=weight.1)
			if(sum(is.na(mu.1.s0))>0){
  				print(paste("Note: ", sum(is.na(mu.1.s0)), " values extrapolated."))
    		c.mat = cbind(szero, mu.1.s0)
    		for(o in 1:length(mu.1.s0)) {
    			if(is.na(mu.1.s0[o])){
    				distance = abs(s0.new - s0.new[o])
    				c.temp = cbind(c.mat, distance)
    				c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where mean is not na
    				new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    				mu.1.s0[o] = new.est[1]   #in case there are multiple matches
    				}
}}
		delta.s = sum((1/weight.0)*mu.1.s0)/n0.all-sum((1/weight.0)*yzero)/n0.all
	}
	return(delta.s)
}
##########################################################
#expit function#
##########################################################
expit = function(x){
	return((exp(x))/(1+exp(x)))
}

##########################################################
#Let's do one iteration#
##########################################################

n1 = 1000
n0=1000

set.seed(1)
data = gen.data(setting=1, n1=n1, n0=n0) 

s1 = data$s1
y1 = data$y1
s0 = data$s0
y0 = data$y0


##########################################################
#observed total treatment effect#
##########################################################
delta.obs = mean(y1)-mean(y0)


##########################################################
#estimate delta.s with observed data - proposed appproach#
##########################################################
delta.s = delta.s.single(sone=s1, szero=s0, yone = y1, yzero = y0)
print(delta.s)
R= 1-delta.s/delta.obs
print(R)
#pretty good

##########################################################
#Let's do 100 iterations#
##########################################################
#only going to look at bias and ese, we can worry about SE and coverage later

set.seed(4)
reps=100
delta.obs.vec = vector(length = reps)   #holds observed delta, simple difference in mean
delta.s.vector =  vector(length = reps)   #holds delta.s 
R.vector= vector(length = reps)   #holds R 

for(i in 1:reps) {
data = gen.data(setting=1, n1=n1, n0=n0)

s1 = data$s1
y1 = data$y1
s0 = data$s0
y0 = data$y0


##########################################################
#observed total treatment effect#
##########################################################
delta.obs = mean(y1)-mean(y0)
delta.obs.vec[i] = delta.obs


##########################################################
#estimate delta.s with observed data - proposed appproach#
##########################################################
delta.s = delta.s.single(sone=s1, szero=s0, yone = y1, yzero = y0)
R= 1-delta.s/delta.obs
delta.s.vector[i] = delta.s
R.vector[i] = R
print(i)
}

#make table
truth = get.true(setting=1)
delta.s.true= truth$delta.s.true
delta.true = truth$delta.true
R.true =  truth$R.true

tab = rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec),mean(delta.s.vector), mean(R.vector)), c(mean(delta.obs.vec)-delta.true,mean(delta.s.vector)-delta.s.true, mean(R.vector)-R.true), c(sd(delta.obs.vec),sd(delta.s.vector), sd(R.vector)))
colnames(tab) = c("Delta","Delta.s","R")
rownames(tab) = c("Truth","Estimate","Bias","ESE")
tab

##########################################################
#Now I am going to add missingness, and do complete case#
##########################################################
#m1 is MCAR
#m2 is MAR
#m3 is MNAR
setting.set = 1

set.seed(4)
reps=100
delta.obs.vec = matrix(nrow = reps, ncol=6)   #holds observed delta, simple difference in mean
delta.s.vector =  matrix(nrow = reps, ncol=6)   #holds delta.s 
R.vector= matrix(nrow = reps, ncol=6)    #holds R 

for(i in 1:reps) {
data = gen.data(setting=setting.set, n1=n1, n0=n0)

s1 = data$s1
y1 = data$y1
s0 = data$s0
y0 = data$y0
m1.1 = rbinom(n1,1,0.70)
m1.0 = rbinom(n0,1,0.70)
#prob = expit(0.02*y1); mean(prob)
#prob = rep(0.7,n1)
#m2.1 = sapply(prob, function(x) rbinom(1,1,x))
#prob = expit(0.03*y0); mean(prob)
#m2.0 = sapply(prob, function(x) rbinom(1,1,x))
#prob = expit(0.15*s1); mean(prob)
#prob = rep(0.7,n1)
#m3.1 = sapply(prob, function(x) rbinom(1,1,x))
#prob = expit(0.12*s0); mean(prob)
#m3.0 = sapply(prob, function(x) rbinom(1,1,x))

m2.1 = rep(1,n1)
m2.1[y1 < median(y1)] = rbinom(sum(y1 < median(y1)),1,0.70)
m2.0 = rep(1,n0)
m2.0[y0 > median(y0)] = rbinom(sum(y0 > median(y0)),1,0.60)
m3.1 = rep(1,n1)
m3.1[s1 < median(s1)] = rbinom(sum(s1 < median(s1)),1,0.70)
m3.0 = rep(1,n0)
m3.0[s0 > median(s0)] = rbinom(sum(s0 > median(s0)),1,0.60)

#known MAR weights
w.y1 = rep(1,n1)
w.y1[y1 < median(y1)] = 0.7
w.y0 = rep(1,n0)
w.y0[y0 > median(y0)] = 0.6
w.s1 = rep(1,n1)
w.s1[s1 < median(s1)] = 0.7
w.s0 = rep(1,n0)
w.s0[s0 > median(s0)] = 0.6

##########################################################
#no missingness#
##########################################################
delta.obs.vec[i,1] = mean(y1)-mean(y0)
delta.s = delta.s.single(sone=s1, szero=s0, yone = y1, yzero = y0)
delta.s.vector[i,1] = delta.s
R.vector[i,1]= 1-delta.s.vector[i,1]/delta.obs.vec[i,1]

##########################################################
#MCAR#
##########################################################
delta.obs.vec[i,2] = mean(y1[m1.1==1])-mean(y0[m1.0==1])
delta.s = delta.s.single(sone=s1[m1.1==1], szero=s0[m1.0==1], yone = y1[m1.1==1], yzero = y0[m1.0==1])
delta.s.vector[i,2] = delta.s
R.vector[i,2]= 1-delta.s.vector[i,2]/delta.obs.vec[i,2]

##########################################################
#MAR#
##########################################################
delta.obs.vec[i,3] = mean(y1[m2.1==1])-mean(y0[m2.0==1])
delta.s = delta.s.single(sone=s1[m2.1==1], szero=s0[m2.0==1], yone = y1[m2.1==1], yzero = y0[m2.0==1])
delta.s.vector[i,3] = delta.s
R.vector[i,3]= 1-delta.s.vector[i,3]/delta.obs.vec[i,3]

##########################################################
#MAR corrected#
##########################################################
delta.obs.vec[i,4] = sum((1/w.y1[m2.1==1])*y1[m2.1==1])/length(y1)-sum((1/w.y0[m2.0==1])*y0[m2.0==1])/length(y0)
delta.s = delta.s.single(sone=s1[m2.1==1], szero=s0[m2.0==1], yone = y1[m2.1==1], yzero = y0[m2.0==1], weight.1 = w.y1[m2.1==1], weight.0 = w.y0[m2.0==1], n0.all = length(y0))
delta.s.vector[i,4] = delta.s
R.vector[i,4]= 1-delta.s.vector[i,4]/delta.obs.vec[i,4]


##########################################################
#MNAR#
##########################################################
delta.obs.vec[i,5] = mean(y1[m3.1==1])-mean(y0[m3.0==1])
delta.s = delta.s.single(sone=s1[m3.1==1], szero=s0[m3.0==1], yone = y1[m3.1==1], yzero = y0[m3.0==1])
delta.s.vector[i,5] = delta.s
R.vector[i,5]= 1-delta.s.vector[i,5]/delta.obs.vec[i,5]

##########################################################
#MNAR weighted#
##########################################################
delta.obs.vec[i,6] = sum((1/w.y1[m3.1==1])*y1[m3.1==1])/length(y1)-sum((1/w.y0[m3.0==1])*y0[m3.0==1])/length(y0)
delta.s = delta.s.single(sone=s1[m3.1==1], szero=s0[m3.0==1], yone = y1[m3.1==1], yzero = y0[m3.0==1], weight.1 = w.y1[m3.1==1], weight.0 = w.y0[m3.0==1], n0.all = length(y0))
delta.s.vector[i,6] = delta.s
R.vector[i,6]= 1-delta.s.vector[i,6]/delta.obs.vec[i,6]



print(i)
}

#make table
truth = get.true(setting=setting.set)
delta.s.true= truth$delta.s.true
delta.true = truth$delta.true
R.true =  truth$R.true

#no missing
tab.no = round(rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec[,1]),mean(delta.s.vector[,1]), mean(R.vector[,1])), c(mean(delta.obs.vec[,1])-delta.true,mean(delta.s.vector[,1])-delta.s.true, mean(R.vector[,1])-R.true), c(sd(delta.obs.vec[,1]),sd(delta.s.vector[,1]), sd(R.vector[,1]))),3)
colnames(tab.no) = c("Delta","Delta.s","R")
rownames(tab.no) = c("Truth","Estimate","Bias","ESE")
#mcar
tab.1 = round(rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec[,2]),mean(delta.s.vector[,2]), mean(R.vector[,2])), c(mean(delta.obs.vec[,2])-delta.true,mean(delta.s.vector[,2])-delta.s.true, mean(R.vector[,2])-R.true), c(sd(delta.obs.vec[,2]),sd(delta.s.vector[,2]), sd(R.vector[,2]))),3)
colnames(tab.1) = c("Delta","Delta.s","R")
rownames(tab.1) = c("Truth","Estimate","Bias","ESE")
#mar
tab.2 = round(rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec[,3]),mean(delta.s.vector[,3]), mean(R.vector[,3])), c(mean(delta.obs.vec[,3])-delta.true,mean(delta.s.vector[,3])-delta.s.true, mean(R.vector[,3])-R.true), c(sd(delta.obs.vec[,3]),sd(delta.s.vector[,3]), sd(R.vector[,3]))),3)
colnames(tab.2) = c("Delta","Delta.s","R")
rownames(tab.2) = c("Truth","Estimate","Bias","ESE")
#mar - fixed
tab.3 = round(rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec[,4]),mean(delta.s.vector[,4]), mean(R.vector[,4])), c(mean(delta.obs.vec[,4])-delta.true,mean(delta.s.vector[,4])-delta.s.true, mean(R.vector[,4])-R.true), c(sd(delta.obs.vec[,4]),sd(delta.s.vector[,4]), sd(R.vector[,4]))),3)
colnames(tab.3) = c("Delta","Delta.s","R")
rownames(tab.3) = c("Truth","Estimate","Bias","ESE")
#mnar
tab.4 = round(rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec[,5]),mean(delta.s.vector[,5]), mean(R.vector[,5])), c(mean(delta.obs.vec[,5])-delta.true,mean(delta.s.vector[,5])-delta.s.true, mean(R.vector[,5])-R.true), c(sd(delta.obs.vec[,5]),sd(delta.s.vector[,5]), sd(R.vector[,5]))),3)
colnames(tab.4) = c("Delta","Delta.s","R")
rownames(tab.4) = c("Truth","Estimate","Bias","ESE")
#mnar - weighted, not fixed
tab.5 = round(rbind(c(delta.true,delta.s.true, R.true), c(mean(delta.obs.vec[,6]),mean(delta.s.vector[,6]), mean(R.vector[,6])), c(mean(delta.obs.vec[,6])-delta.true,mean(delta.s.vector[,6])-delta.s.true, mean(R.vector[,6])-R.true), c(sd(delta.obs.vec[,6]),sd(delta.s.vector[,6]), sd(R.vector[,6]))),3)
colnames(tab.5) = c("Delta","Delta.s","R")
rownames(tab.5) = c("Truth","Estimate","Bias","ESE")
big.table = rbind(c("No missing", "",""),tab.no, c("MCAR", "",""), tab.1, c("MAR","",""), tab.2, c("MAR-weighted", "",""), tab.3, c("MNAR", "",""), tab.4, c("MNAR - weighted", "",""), tab.5)
big.table
library(quantreg)
latex.table(big.table, "sim_missingness", row.names = T, col.names = F, caption = "", dcolumn = T)


delta.obs.vec=as.data.frame(delta.obs.vec)
delta.s.vector = as.data.frame(delta.s.vector)
R.vector = as.data.frame(R.vector)
names(delta.obs.vec) = names(delta.s.vector) = names(R.vector) = c("No missing", "MCAR", "MAR","MAR-w","MNAR", "MNAR-w")
par(mfrow=c(1,3))
boxplot(delta.obs.vec, main="Delta")
abline(h=delta.true, col="red", lwd=2, lty=2) 
boxplot(delta.s.vector, main = "Delta.s")
abline(h=delta.s.true, col="red", lwd=2, lty=2) 
boxplot(R.vector, main = "R")
abline(h=R.true, col="red", lwd=2, lty=2) 

##########################################################
#Some pictures of the conditional mean thinking about extrapolating#
##########################################################

#### Let's look at this conditional mean function 
#make some data
set.seed(1)
setting.set=2
data = gen.data(setting=setting.set, n1=n1, n0=n0)

s1 = data$s1
y1 = data$y1
s0 = data$s0
y0 = data$y0

#just working with the treated group right now
#this is the s "grid" you care about, say s=s1
#plot1
s.grid=s1
mu.1.s0 = sapply(s.grid,pred.smooth,zz=s1, y1=y1)
plot(s.grid, y1, col="light pink",pch=20) #this is the true data only because s.grid=s1 
points(s.grid,mu.1.s0, pch=20)
legend("bottomright", pch=20, col = c("black","light pink"), c("Kernel estimate","True data"))
###plot 2 
s.grid=s1
mu.1.s0 = sapply(s.grid,pred.smooth,zz=s1, y1=y1)
plot(s.grid, y1, col="light pink",pch=20) #this is the true data only because s.grid=s1 
points(s.grid,mu.1.s0, pch=20)
#now let's color the ones that are missing as blue
m2.1 = rep(1,n1)
m2.1[y1>40] = 1-rbinom(sum(y1>40),1,0.80)
points(s.grid[m2.1==0],y1[m2.1==0], pch=20, col = "light blue")
legend("bottomright", pch=20, col = c("black","light pink","light blue"), c("Kernel estimate","True data", "Mising"))

#plot 3, I guess we would sort of be flipping this and trying to predict s from y, 
mu.1.s0 = sapply(y1,pred.smooth,zz=y1[m2.1==1], y1=s1[m2.1==1])
plot(y1, s1, col="light pink",pch=20) #this is the true data only because s.grid=s1 
points(y1,mu.1.s0, pch=20)
points(y1[m2.1==0],s1[m2.1==0], pch=20, col = "light blue")

legend("bottomright", pch=20, col = c("black","light pink"), c("Kernel estimate","True data"))
