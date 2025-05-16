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