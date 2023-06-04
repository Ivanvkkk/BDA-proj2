# question 1
require(insuranceData)
data(dataCar)
#You may need to set the working directory first before loading the dataset
#setwd("Desktop/S2/Bayesian Data Analysis/assignment-2")
#The first 6 rows of the dataframe
print.data.frame(dataCar[1:6,])

library(INLA)
str(dataCar)
# Center and scale the non-categorical covariates
dataCar$veh_value <- scale(dataCar$veh_value)[,1]
dataCar$exposure <- scale(dataCar$exposure)[,1]
# Convert integers to factors for categorical covariates
dataCar$veh_age <- as.factor(dataCar$veh_age)
dataCar$agecat <- as.factor(dataCar$agecat)

str(dataCar)

#1.a
#评判标准：log marginal likelihood, DIC and NLSCPO 
#model$mlik[1]
#model$cpo$cpo
#model$dic$dic

formula_clm <- clm ~ veh_value + exposure + veh_body + veh_age + gender + area + agecat
#priors.beta <- list(mean.intercept = 0, prec.intercept = 100,
#               mean = list('veh_value' = 0, 'exposure' = 0),
#               prec = list('veh_value' = 100, 'exposure' = 100))

prior.beta <- list(mean.intercept = 0, prec.intercept = 0.001,
                   mean = 0, prec = 0.001)

model <- inla(formula_clm, family = "binomial", data = dataCar, 
              control.fixed = prior.beta,
              control.family = list(link = "logit"),
              control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
              )
summary(model)
model$summary.fitted.values$mean


# 1.b
formula_numclaims <- numclaims ~ veh_value + exposure + veh_body + veh_age + gender + area + agecat

prior.beta <- list(mean.intercept = 0, prec.intercept = 1/log(4)^2,
                   mean = 0, prec = 1/(log(13)/6)^2)
model_poisson <- inla(formula_numclaims, family = "poisson", data = dataCar, 
                      control.fixed = prior.beta, 
                      control.family = list(link = "log"),
                      control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
                      )
summary(model_poisson)

# 1.c
formula_numclaims <- numclaims ~ veh_value + exposure + veh_body + veh_age + gender + area + agecat
prior.beta <- list(mean.intercept = 0, prec.intercept = 1/log(4)^2,
                   mean = 0, prec = 1/(log(13)/6)^2)
model_zero <- inla(formula_numclaims, 
                        family = "zeroinflatedpoisson1", 
                        data = dataCar, 
                        control.fixed = prior.beta,
                        control.family = list(link = "log"),
                        control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)
summary(model_zero)

# 1.d
#lec5 ,p62
# sigma_obs:
sigma.unif.prior = "expression:
b = 20;
log_dens= (theta>=(-2*log(b)))*(-log(b)-theta/2-log(2))+
(theta<(-2*log(b)))*(-Inf); return(log_dens);"
#sigma_alpha:
sigma.unif.prior.random.eff = "expression:
b = 20;
log_dens = (theta>=(-2*log(b)))*(-log(b)-theta/2-log(2)) +
(theta<(-2*log(b)))*(-Inf); return(log_dens);"
b=20;
prec.prior <- list(prec=list(prior = sigma.unif.prior,
                             initial = -2*log(b)+1,fixed = FALSE))
prec.prior.random.eff <- list(prec=list(prior =
                                          sigma.unif.prior.random.eff, 
                                        initial = -2*log(b)+1, fixed = FALSE))
formula_random <- numclaims ~ 
  veh_value + exposure + veh_age + gender + area + agecat +
  veh_value:exposure + veh_age:gender + area:gender + 
  f(veh_body, model = "iid",hyper= prec.prior.random.eff) 

prior.beta <- list(mean.intercept = 0, prec.intercept = 1/log(4)^2,
                   mean = 0, prec = 1/(log(13)/6)^2)

model_random <- inla(formula_random, family = "poisson", data = dataCar, 
                     control.fixed=prior.beta,
                     control.inla = list(strategy = "laplace", npoints = 40),
                     control.compute = list(config = TRUE,dic = TRUE, cpo=TRUE))
summary(model_random)


# 1.e
nbsamp=1000
n=nrow(dataCar)
yrep1 = matrix(0,nrow=n,ncol=nbsamp)
yrep2 = matrix(0,nrow=n,ncol=nbsamp)
yrep3 = matrix(0,nrow=n,ncol=nbsamp)

poisson.samples=inla.posterior.sample(nbsamp, result=model_poisson)
zero.samples=inla.posterior.sample(nbsamp, result=model_zero)
random.samples=inla.posterior.sample(nbsamp, result=model_random)
predictor.samples.poisson=inla.posterior.sample.eval(function(...) {Predictor},
                                                     poisson.samples)
predictor.samples.zero=inla.posterior.sample.eval(function(...) {Predictor},
                                                  zero.samples)
predictor.samples.random=inla.posterior.sample.eval(function(...) {Predictor},
                                                    random.samples)

for (row.num in 1:n){   
  yrep1[row.num,]<- rpois(nbsamp,
                           lambda=exp(predictor.samples.poisson[row.num,]))
  yrep2[row.num,]<- rpois(nbsamp,
                           lambda=exp(predictor.samples.zero[row.num,]))
  yrep3[row.num,]<- rpois(nbsamp,
                           lambda=exp(predictor.samples.random[row.num,]))
}

plot.post.pred.test <- function(yrep) {
  par(mfrow = c(3, 2))
  par(mar = c(1.7, 1.7, 1.7, 1.7))
  
  for (i in 0:4) {
    numclaims <- apply(yrep, 2, function(x) sum(x == i))
    hist(numclaims, xlim = c(min(numclaims), max(numclaims)), 
         main = paste("Numclaims =", i), xlab = "Value", ylab = "Frequency")
    abline(v = sum(dataCar$numclaims == i), col = 'red', lwd = 2)
  }
  
  par(mfrow = c(1, 1))
}

plot.post.pred.test(yrep1)
plot.post.pred.test(yrep2)
plot.post.pred.test(yrep3)

y_hat_poisson = model_poisson$summary.fitted.values[,1]
y_hat_zero = model_zero$summary.fitted.values[,1]
y_hat_random = model_random$summary.fitted.values[,1]
rmse_poisson <- sqrt(mean((y_hat_poisson - dataCar$numclaims)^2))
rmse_zero <- sqrt(mean((y_hat_zero - dataCar$numclaims)^2))
rmse_random <- sqrt(mean((y_hat_random - dataCar$numclaims)^2))


# Question2
#2.a
study<-read.csv("Barcelona.csv")
head(study)
#remove the missing value
study <- study[apply(study != "", 1, all), ]
#calculate the mean and sd
vars_to_scale <- c('no2gps_24h', 'maxwindspeed_24h', 'precip_24h', 'sec_noise55_day', 'age_yrs', 'tmean_24h')
mean_sd_list <- lapply(vars_to_scale, function(var) {
  list(mean = mean(study[[var]]), sd = sd(study[[var]]))
})
# Center and scale the non-categorical covariates
study$no2gps_24h <- scale(study$no2gps_24h)[,1]
study$maxwindspeed_24h <- scale(study$maxwindspeed_24h)[,1]
study$precip_24h <- scale(study$precip_24h)[,1]
study$sec_noise55_day <- scale(study$sec_noise55_day)[,1]
study$age_yrs <- scale(study$age_yrs)[,1]
study$tmean_24h <- scale(study$tmean_24h)[,1]

# Convert integers to factors for categorical covariates
cols_to_factor <- c("gender", "on_a_diet", "alcohol", "drugs", "sick", "other_factors", 
                    "education", "smoke", "access_greenbluespaces_300mbuff")
study[cols_to_factor] <- lapply(study[cols_to_factor], as.factor)

str(study)

formula<- log(stroop_test_performance) ~ gender+on_a_diet+alcohol+drugs+sick+other_factors+education+
  smoke+no2gps_24h+maxwindspeed_24h+precip_24h+sec_noise55_day+
  access_greenbluespaces_300mbuff+age_yrs+tmean_24h
model_linear <- inla(formula, family = "gaussian", data =study, 
                      control.fixed = list(mean = 0, prec = 1/100), 
                      control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)
summary(model_linear)

#2.b
formula_poisson<- sadness ~ gender+on_a_diet+alcohol+drugs+sick+other_factors+education+
  smoke+no2gps_24h+maxwindspeed_24h+precip_24h+sec_noise55_day+
  access_greenbluespaces_300mbuff+age_yrs+tmean_24h

model_poisson <- inla(formula_poisson, family = "poisson", data =study, 
                      control.fixed = list(mean = 0, prec = 1/100), 
                      control.family = list(link = "log"),
                      control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)
summary(model_poisson)

#2.c
study$Person_ID <- as.factor(study$Person_ID)
formula.linear.random<- log(stroop_test_performance) ~ 
  gender+on_a_diet+alcohol+drugs+sick+other_factors+education+
  smoke+no2gps_24h+maxwindspeed_24h+precip_24h+sec_noise55_day+
  access_greenbluespaces_300mbuff+age_yrs+tmean_24h+
  f(Person_ID, model = "iid") 
formula.poisson.random<- sadness ~ gender+on_a_diet+alcohol+
  drugs+sick+other_factors+education+
  smoke+no2gps_24h+maxwindspeed_24h+precip_24h+sec_noise55_day+
  access_greenbluespaces_300mbuff+age_yrs+tmean_24h+
  f(Person_ID, model = "iid") 

model.linear.random <- inla(formula.linear.random, family = "gaussian", data =study, 
                     control.fixed = list(mean = 0, prec = 1/100), 
                     control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)

model.poisson.random <- inla(formula.poisson.random, family = "poisson", data =study, 
                      control.fixed = list(mean = 0, prec = 1/100), 
                      control.family = list(link = "log"),
                      control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)

summary(model.linear.random)
summary(model.poisson.random)
#2.d
nbsamp=1000
n=nrow(study)
yrep1 = matrix(0,nrow=n,ncol=nbsamp)
yrep2 = matrix(0,nrow=n,ncol=nbsamp)

poisson.samples=inla.posterior.sample(n=nbsamp, result=model_poisson)
random.samples=inla.posterior.sample(n=nbsamp, result=model.poisson.random)

predictor.samples.poisson=inla.posterior.sample.eval(function(...) {Predictor},
                                                     poisson.samples)
predictor.samples.random=inla.posterior.sample.eval(function(...) {Predictor},
                                                    random.samples)

for (row.num in 1:n){   
  yrep1[row.num,]<- rpois(n=nbsamp,
                          lambda=exp(predictor.samples.poisson[row.num,]))
  yrep2[row.num,]<- rpois(n=nbsamp,
                          lambda=exp(predictor.samples.random[row.num,]))
}

plot.post.pred.test<-function(yrep){
  sadness.per.gender.samples=aggregate(yrep,list(study$gender), mean)
  sadness.per.gender.in.data=aggregate(study$sadness,list(study$gender), mean)
  par(mfrow=c(3,1))
  par(mar=c(1.7,1.7,1.7,1.7))
  for(it in 1:3) { 
    x=as.numeric(sadness.per.gender.samples[it,2:(nbsamp+1)])
    sadness.on.data=sadness.per.gender.in.data[it,2]
    xmin=min(min(x),sadness.on.data)
    xmax=max(max(x),sadness.on.data)
    hist(as.numeric(sadness.per.gender.samples[it,2:(nbsamp+1)]),
       col="gray40",main=sadness.per.gender.samples[it,1],xlim=c(xmin,xmax))
    abline(v=sadness.per.gender.in.data[it,2],col="red",lwd=2) 
  }
  par(mfrow=c(1,1))
}
plot.post.pred.test(yrep1)
plot.post.pred.test(yrep2)

y_hat_poisson = model_poisson$summary.fitted.values[,1]
y_hat_random = model.poisson.random$summary.fitted.values[,1]
rmse_poisson <- sqrt(mean((y_hat_poisson - study$sadness)^2))
rmse_random <- sqrt(mean((y_hat_random - study$sadness)^2))

#2.e
# create the new data
new_person <- data.frame(
  Person_ID = factor("286"),
  gender = factor("Woman"),
  on_a_diet = factor("Yes"),
  alcohol = factor("No"),
  drugs = factor("No"),
  sick = factor("No"),
  other_factors = factor("No"),
  education = factor("University"),
  smoke = factor("Yes"),
  no2gps_24h = 80,
  maxwindspeed_24h = 10,
  precip_24h = 50,
  sec_noise55_day = 10000,
  access_greenbluespaces_300mbuff = factor("Yes"),
  age_yrs = 40,
  tmean_24h = 25,
  stroop_test_performance = NA,
  sadness = NA
)

# select the used variable
select.var = c('gender', 'on_a_diet', 'alcohol', 'drugs', 'sick', 'other_factors',
'education', 'smoke', 'no2gps_24h', 'maxwindspeed_24h', 'precip_24h', 'sec_noise55_day',
'access_greenbluespaces_300mbuff', 'age_yrs', 'tmean_24h','Person_ID','stroop_test_performance','sadness')
study.sub = study[, select.var]

# scale the new data
mean_array <- c()
sd_array <- c()

for (i in 1:length(mean_sd_list)) {
  mean_array[i] <- mean_sd_list[[i]]$mean
  sd_array[i] <- mean_sd_list[[i]]$sd
}
mean_sd_df <- data.frame(variable = vars_to_scale, Mean = mean_array, SD =sd_array )
for (i in 1:nrow(mean_sd_df)) {
  var <- mean_sd_df$variable[i]
  new_person[[var]] <- (new_person[[var]] - mean_sd_df$Mean[i]) / mean_sd_df$SD[i]
}
#combine the data
study.sub = rbind(study.sub, new_person)

# use the model with new data
model.linear.random.new <- inla(formula.linear.random, family = "gaussian", data =study.sub, 
                            control.fixed = list(mean = 0, prec = 1/100), 
                            control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)

model.poisson.random.new <- inla(formula.poisson.random, family = "poisson", data =study.sub, 
                             control.fixed = list(mean = 0, prec = 1/100), 
                             control.family = list(link = "log"),
                             control.compute = list(config = TRUE,cpo=TRUE, dic = TRUE)
)
study.linear.samp=inla.posterior.sample(n=nbsamp, result=model.linear.random.new ,selection= list(Predictor=nrow(study.sub)))
study.poisson.samp=inla.posterior.sample(n=nbsamp, result=model.poisson.random.new,selection= list(Predictor=nrow(study.sub)))

predictor.linear.samples=exp(unlist(lapply(study.linear.samp, function(x)(x$latent[1]))))
predictor.poisson.samples=exp(unlist(lapply(study.poisson.samp, function(x)(x$latent[1]))))
library(ggplot2)

# Plot estimated density for stroop_test_performance
ggplot() +
  geom_density(aes(x = predictor.linear.samples), color = "blue", fill = "blue", alpha = 0.5) +
  labs(title = "Posterior Predictive Distribution for Stroop Test Performance",
       x = "Stroop Test Performance",
       y = "Density")

# Plot histogram for sadness
ggplot() +
  geom_histogram(aes(x = predictor.poisson.samples), color = "red", fill = "red", alpha = 0.5, bins = 15) +
  labs(title = "Posterior Predictive Distribution for Sadness",
       x = "Sadness",
       y = "Frequency")

