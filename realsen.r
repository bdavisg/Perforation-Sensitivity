library(sensitivity)
library(rstan)

setwd("C:/Users/brad.davis/Desktop/Dakota/Sens")

#General Stan Settings
rstan_options(auto_write = TRUE) #saves compiled model as an .rds automatically
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

#load processed data and previous Stan file
source('sim.data.r')
source('post.data.r')
fit_predict <- readRDS("fit_predict2.rds") #this will load the fit model that I already ran
samples<-extract(fit_predict)
rho<-mean(samples$rho)
alpha<-mean(samples$alpha)
sigma<-mean(samples$sigma)

#Create vector x1 and x2 from sim.data.r
X1<-stan_data$x1[1:145,]
X2<-stan_data$x1[146:290,]

#Create vector X of samples to simulation from fixed GP model. Called x2 in stan model.
sobol<-sobol2002(model=NULL,X1,X2,nboot=1000)
x2<-sobol$X

#Create new vector of data for Stan model to run
N2<-nrow(x2)
D2<-ncol(x2)

#stan_data1 <- list(N1=stan_data$N1,D1=stan_data$D1,N2=N2,D2=D2,x1=stan_data$x1,y1=stan_data$y1,x2=x2,rho=rho,alpha=alpha,sigma=sigma)
stan_data1 <- list(N1=stan_data$N1,D1=stan_data$D1,N2=stan_data$N2,D2=stan_data$D2,x1=stan_data$x1,y1=stan_data$y1,x2=stan_data$x2,rho=rho,alpha=alpha,sigma=sigma)

#Run stan model to generate samples
fit_sens <- stan(file="perf-sens.stan",data=stan_data1,iter=1000, seed=249, refresh=10, control = list(max_treedepth=10))

#Complete Sensitivity Analysis
y2 <- extract(fit_sens)$y2[1,] #extract x2 results from model
tell(sobol,y2)
print(x1)
plot(x1)