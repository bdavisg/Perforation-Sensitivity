# This is a script to use a fixed gaussian process emulator in a global Sobol 
# sensitivity analysis using the Sensitivity package

#Load required libraries
library(sensitivity)
library(rethinking)

#Set working directory
setwd("C:/Users/brad.davis/Desktop/Dakota/Sens")

#Load Utilities from Betancourt for post-processing
util <- new.env() #Not sure that I need this line. Looks like it creates a clean workspace
source('stan_utility.R', local=util)

#General Stan Settings
rstan_options(auto_write = TRUE) #saves compiled model as an .rds automatically
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

#load processed data and previous Stan file
source('sim.data.r')
fit_predict <- readRDS("fit_predict2.rds") #this will load the fit model that I already ran

#Extract samples and calculate the mean for GP model parameters to use in fixed GP model
samples<-extract(fit_predict)
rho<-mean(samples$rho)
alpha<-mean(samples$alpha)
sigma<-mean(samples$sigma)

#Create vector x1 and x2 from sim.data.r
#Must use at least 226 data points or Senitivity package will not work
X1<-as.data.frame(stan_data$x1[1:226,])
X2<-as.data.frame(stan_data$x1[226:1,])

#Create vector X of samples to simulation from fixed GP model. Called x2 in stan model.
sens<-sobol2007(model=NULL,X1,X2,nboot=1000)
x2<-sens$X #Becomes vector x2 in Stan model

#Create new vector of data for Stan model to run
N2<-nrow(x2)
D2<-ncol(x2)

#Load data into a single list for Stan
stan_data1 <- list(N1=stan_data$N1,D1=stan_data$D1,N2=N2,D2=D2,x1=stan_data$x1,y1=stan_data$y1,x2=x2,rho=rho,alpha=alpha,sigma=sigma)
#stan_data1 <- list(N1=stan_data$N1,D1=stan_data$D1,N2=stan_data$N2,D2=stan_data$D2,x1=stan_data$x1,y1=stan_data$y1,x2=stan_data$x2,rho=rho,alpha=alpha,sigma=sigma)

#Run stan model to generate samples
fit_sens <- stan(file="perf-sens.stan",data=stan_data1,iter=2000, seed=249, refresh=10, control = list(max_treedepth=11))

#Check model performance
util$check_all_diagnostics(fit_sens)

#Save Stan model results for use in future sessions
saveRDS(fit_sens,"fit_sens1.rds")

#Extract results for y2
y2 <- extract(fit_sens)$y2[1,]

#Plot results to make sure that data looks reasonable
plot(stan_data$x1,stan_data$y1)
points(stan_data$x2,y2,col='black',pch=16, cex=1)

#Complete Sensitivity Analysis
tell(sens,y2)
print(sens)
plot(sens)
