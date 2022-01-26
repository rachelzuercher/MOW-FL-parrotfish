# Load libraries ####
library(dismo)
library(gbm)
library(reshape2)
library(plyr)
library(dplyr)
library(xlsx)
library(ggplot2)
library(stringr)
library(ggBRT)
library(ape)
library(scales)
library(data.table)
library(fmsb)

options(scipen=999)

# Install ggBRT 
#install.packages("devtools") # if package "devtools" not already installed
#devtools::install_github("JBjouffray/ggBRT")
#library(ggBRT)

# Set working directory and load data ####
setwd("")
rvc_pf <- read.csv("RVC_parrotfish_2012-2018.csv", sep=",", quote="", na.strings=c("NA", "", " "), header=TRUE)


# Set up new dataframe to deal with spatial autocorrelation in Midnight and Blue presence/absence models ####

# take every third site by latitude
rvc_pf3 <- rvc_pf[order(rvc_pf$Latitude),]
rvc_pf3 <- rvc_pf3[-seq(3, NROW(rvc_pf3), by = 3),]



# Test for collinearity among variables and calculation of VIFs ####

# Select the set of predictors from the rvc_pf dataframe
predictors <- c(12:44)
predictors.numeric <- c(12,14:21,23:26,28:36,37:44)

# Assess collinearity with Pearson correlation coefficient
library(corrplot)
coeff <- cor(rvc_pf[,predictors.numeric], method="pearson",use="pairwise.complete.obs")
coeff[which(abs(coeff[]) > 0.8)] # no correlation > 0.8

# Correlation plot
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(coeff,method="number", col=col(200),tl.cex=0.8,number.cex = 0.7,
         tl.col="black", tl.srt=45)
# to save, use the 'Export as PDF' option in the Plots window 18x15 

predictors.posttrim <- c(12,16:21,23:25,28,31:36,37:44)
coeff1 <- cor(rvc_pf[,predictors.posttrim], method="pearson",use="pairwise.complete.obs")
coeff1[which(abs(coeff1[]) > 0.8)] # no correlation > 0.8
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot(coeff1,method="number", col=col(200),tl.cex=0.8,number.cex = 0.7,
         tl.col="black", tl.srt=45)
# to save, use the 'Export as PDF' option in the Plots window 18x15 

# Assess collinearity with Variable Inflation Factor (VIF)
#source("vif_func.R") # available at https://gist.github.com/fawda123/4717702#file-vif_fun-r

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  library(fmsb)
  
  if(any(!'data.frame' %in% class(in_frame))) in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    return(names(in_dat))
  }
}

vif_func(in_frame=rvc_pf[,predictors.posttrim],thresh=8,trace=T) 



# Make predictor lists for each species ####

pred.midnightPA <- c(12,13,16:25,27:28,31:37,38,42:44)
pred.midnight <- c(12,13,16:25,27:28,31:37,38,42:44)

pred.bluePA <- c(12,13,16:25,27:28,31:37,39,42:44)
pred.blue <- c(12,13,16:25,27:28,31:37,39,42:44)

pred.rainbowPA <- c(12,13,16:25,27:28,31:37,40,42:44)
pred.rainbow <- c(12,13,16:25,27:28,31:37,40,42:44)

pred.stoplightPA <- c(12,13,16:25,27:28,31:37,41,42:44)
pred.stoplight <- c(12,13,16:25,27:28,31:37,41,42:44)



# Sample sizes ####

sample.sizes <- data.frame(Species=c("Midnight","Blue","Rainbow","Stoplight"),
                           n_Present=c(count(rvc_pf[rvc_pf$midnightPA==1,])[[1]],
                                       count(rvc_pf[rvc_pf$bluePA==1,])[[1]],
                                       count(rvc_pf[rvc_pf$rainbowPA==1,])[[1]], 
                                       count(rvc_pf[rvc_pf$stoplightPA==1,])[[1]]))


pf.sites <- rvc_pf[,c(1,8:11)]
pf.sites$sum <- rowSums(pf.sites[,c(2:5)])
count(pf.sites[pf.sites$sum>0,])



# Test BRT parameters ####

tree_complexity <- 1:5
learning_rate <- c(0.01,0.05,0.001,0.0001)
bag_fraction <- c(0.5,0.75,0.9)

names <- c("tc", "lr", "bf", "tot_dev", "resid_dev", "corr",
           "AUC", "perc_expl", "cv_dev", "cv_corr", "cv_AUC", "cv_perc_expl")

# make a place to put results from the loop below
parameter_tests <- data.frame(matrix(NA, nrow = 1, ncol = 12))
colnames(parameter_tests) <- names

output <- data.frame(matrix(NA, nrow = 1, ncol = 12))
colnames(output) <- names

for (i in 1:length(tree_complexity)) {
  for (j in 1:length(learning_rate)) {
    for (k in 1:length(bag_fraction)) {
      
      model <- gbm.step(data = rvc_pf, 
                        gbm.x = pred.midnightPA, # run this loop for each pres/abs and each biomass model
                        gbm.y = 8,
                        family = "gaussian", 
                        tree.complexity = tree_complexity[i],
                        learning.rate = learning_rate[j],
                        bag.fraction = bag_fraction[k])
      
      output$tc <- tree_complexity[i]
      output$lr <- learning_rate[j]
      output$bf <- bag_fraction[k]
      output$tot_dev <- model$self.statistics$mean.null
      output$resid_dev <- model$self.statistics$mean.resid
      output$corr <- model$self.statistics$correlation
      output$AUC <- model$self.statistics$discrimination
      output$perc_expl <- (1 - model$self.statistics$mean.resid / model$self.statistics$mean.null)*100
      output$cv_dev <- model$cv.statistics$deviance.mean
      output$cv_corr <- model$cv.statistics$correlation.mean
      output$cv_AUC <- model$cv.statistics$discrimination.mean
      output$cv_perc_expl <- (1 - model$cv.statistics$deviance.mean / model$self.statistics$mean.null)*100
      
      parameter_tests <- rbind(parameter_tests, output)
      
    }
  }
}





# BRT models: midnight (Scarus coelestinus) PA ####

set.seed(10)
brt_PA_midnight <- gbm.step(data = rvc_pf, 
                            gbm.x = pred.midnightPA, 
                            gbm.y = 8, 
                            family = "bernoulli", 
                            tree.complexity = 5,
                            learning.rate = 0.01, 
                            bag.fraction = 0.75) 

Perf <- ggPerformance(midnight.PA=brt_PA_midnight) 
round(Perf,2)
ggInfluence(brt_PA_midnight, col.bar = "midnightblue", 
            show.signif=FALSE,
            main="Presence/Absence Midnight PF")

# PA Midnight signif
set.seed(12)
brt_PA_midnight.signif <- gbm.step(data = rvc_pf,
                                   gbm.x = c(16,17,18,22,23,24,33,37,38), 
                                   gbm.y = 8,
                                   family = "bernoulli", 
                                   tree.complexity = 5,
                                   learning.rate = 0.01, 
                                   bag.fraction = 0.75)
Perf <- ggPerformance(midnight.PA=brt_PA_midnight.signif) 
round(Perf,2)
ggInfluence(brt_PA_midnight.signif,col.bar = "midnightblue", 
            show.signif=FALSE,
            main="Presence/Absence Midnight PF")
gbm.plot(brt_PA_midnight.signif, n.plots=9, write.title = FALSE)

brt_PA_midnight.signif.prerun<- plot.gbm.4list(brt_PA_midnight.signif)
brt_PA_midnight.signif.boot <- gbm.bootstrap.functions(brt_PA_midnight.signif, list.predictors=brt_PA_midnight.signif.prerun, n.reps=100)

# Moran's I for PA midnight model #
distancesPA <- as.matrix(dist(cbind(rvc_pf$Latitude, rvc_pf$Longitude)))
distancesPA[1:5, 1:5]
distancesPA.inverse <- 1/distancesPA
diag(distancesPA.inverse) <- 0
distancesPA.inverse[is.infinite(distancesPA.inverse)] <- 0 
distancesPA.inverse[1:5, 1:5]
Moran.I(brt_PA_midnight.signif$residuals, distancesPA.inverse, na.rm=TRUE)


# EVERY 3rd SITE MIDNIGHT P/A
set.seed(12)
brt_PA_midnight.every3 <- gbm.step(data = rvc_pf3, 
                                   gbm.x = pred.midnightPA, 
                                   gbm.y = 8, 
                                   family = "bernoulli", 
                                   tree.complexity = 5,
                                   learning.rate = 0.01, 
                                   bag.fraction = 0.75) 

Perf <- ggPerformance(midnight.PA=brt_PA_midnight.every3) 
round(Perf,2)
ggInfluence(brt_PA_midnight.every3, col.bar = "midnightblue", 
            show.signif=FALSE,
            main="Presence/Absence Midnight PF (Every3)")

set.seed(12)
brt_PA_midnight.every3.signif <- gbm.step(data = rvc_pf3, 
                                          gbm.x = c(16,18,24,33,35,37,38),
                                          gbm.y = 8,
                                          family = "bernoulli", 
                                          tree.complexity = 5,
                                          learning.rate = 0.01, 
                                          bag.fraction = 0.75)
Perf <- ggPerformance(midnight.PA=brt_PA_midnight.every3.signif) 
round(Perf,2)
ggInfluence(brt_PA_midnight.every3.signif,col.bar = "midnightblue", 
            show.signif=FALSE,
            main="Presence/Absence Midnight PF (Every3)")

brt_PA_midnight.every3.signif.prerun<- plot.gbm.4list(brt_PA_midnight.every3.signif)
brt_PA_midnight.every3.signif.boot <- gbm.bootstrap.functions(brt_PA_midnight.every3.signif, list.predictors=brt_PA_midnight.every3.signif.prerun, n.reps=100)

# Moran's I for PA midnight model every 3rd site #
distancesPA <- as.matrix(dist(cbind(rvc_pf3$Latitude, rvc_pf3$Longitude)))
distancesPA[1:5, 1:5]
distancesPA.inverse <- 1/distancesPA
diag(distancesPA.inverse) <- 0
distancesPA.inverse[is.infinite(distancesPA.inverse)] <- 0 
distancesPA.inverse[1:5, 1:5]
Moran.I(brt_PA_midnight.every3.signif$residuals, distancesPA.inverse, na.rm=TRUE)




# BRT models: midnight (Scarus coelestinus) Biomass ####
set.seed(4)
brt_bio_midnight <- gbm.step(data = rvc_pf[rvc_pf$logp1_midnightPF_biomass_kg>0,], 
                             gbm.x = pred.midnight, 
                             gbm.y = 4, 
                             family = "gaussian", 
                             tree.complexity = 5,
                             learning.rate = 0.001, 
                             bag.fraction = 0.75)

Perf <- ggPerformance(midnight.biomass=brt_bio_midnight) 
round(Perf,2)
ggInfluence(brt_bio_midnight, col.bar = "midnightblue", 
            show.signif=FALSE,
            main="Biomass Midnight PF")

# PA Midnight signif
set.seed(4)
brt_bio_midnight.signif <- gbm.step(data = rvc_pf[rvc_pf$logp1_midnightPF_biomass_kg>0,], 
                                    gbm.x = c(13,37,38), 
                                    gbm.y = 4, 
                                    family = "gaussian", 
                                    tree.complexity = 5,
                                    learning.rate = 0.001, 
                                    bag.fraction = 0.75)
Perf <- ggPerformance(midnight.biomass=brt_bio_midnight.signif) 
round(Perf,2)
ggInfluence(brt_bio_midnight.signif,col.bar = "midnightblue", 
            show.signif=FALSE,
            main="Biomass Midnight PF")

brt_bio_midnight.signif.prerun<- plot.gbm.4list(brt_bio_midnight.signif)
brt_bio_midnight.signif.boot <- gbm.bootstrap.functions(brt_bio_midnight.signif, list.predictors=brt_bio_midnight.signif.prerun, n.reps=100)

# Moran's I for Biomass midnight model
distances <- as.matrix(dist(cbind(rvc_pf[rvc_pf$logp1_midnightPF_biomass_kg>0,]$Latitude, 
                                  rvc_pf[rvc_pf$logp1_midnightPF_biomass_kg>0,]$Longitude)))
distances[1:5, 1:5]
distances.inverse <- 1/distances
diag(distances.inverse) <- 0
distances.inverse[1:5, 1:5]
Moran.I(brt_bio_midnight.signif$residuals, distances.inverse, na.rm=TRUE)





# BRT models: blue (Scarus coeruleus) PA ####

set.seed(7)
brt_PA_blue <- gbm.step(data = rvc_pf, 
                        gbm.x = pred.bluePA, 
                        gbm.y = 9, 
                        family = "bernoulli", 
                        tree.complexity = 5,
                        learning.rate = 0.01, 
                        bag.fraction = 0.75) 

Perf <- ggPerformance(blue.PA=brt_PA_blue) 
round(Perf,2)
ggInfluence(brt_PA_blue, col.bar = "#1380BE", 
            show.signif=FALSE,
            main="Presence/Absence Blue PF")

# PA blue signif
set.seed(2)
brt_PA_blue.signif <- gbm.step(data = rvc_pf, 
                               gbm.x = c(13,16,17,18,24,25,33,37,39), 
                               gbm.y = 9, 
                               family = "bernoulli", 
                               tree.complexity = 5,
                               learning.rate = 0.01, 
                               bag.fraction = 0.75)
Perf <- ggPerformance(blue.PA=brt_PA_blue.signif) 
round(Perf,2)
ggInfluence(brt_PA_blue.signif,col.bar = "#1380BE", 
            show.signif=FALSE,
            main="Presence/Absence Blue PF")

brt_PA_blue.signif.prerun<- plot.gbm.4list(brt_PA_blue.signif)
brt_PA_blue.signif.boot <- gbm.bootstrap.functions(brt_PA_blue.signif, list.predictors=brt_PA_blue.signif.prerun, n.reps=100)

# Moran's I for Pres/Abs blue model
distancesPA <- as.matrix(dist(cbind(rvc_pf$Latitude, rvc_pf$Longitude)))
distancesPA[1:5, 1:5]
distancesPA.inverse <- 1/distancesPA
diag(distancesPA.inverse) <- 0
distancesPA.inverse[is.infinite(distancesPA.inverse)] <- 0 
distancesPA.inverse[1:5, 1:5]
Moran.I(brt_PA_blue.signif$residuals, distancesPA.inverse, na.rm=TRUE)


# EVERY 3rd SITE BLUE P/A
set.seed(7)
brt_PA_blue.every3 <- gbm.step(data = rvc_pf3, 
                               gbm.x = pred.bluePA, 
                               gbm.y = 9, 
                               family = "bernoulli", 
                               tree.complexity = 5,
                               learning.rate = 0.01, 
                               bag.fraction = 0.75) 

Perf <- ggPerformance(blue.PA=brt_PA_blue.every3) 
round(Perf,2)
ggInfluence(brt_PA_blue.every3, col.bar = "#1380BE", 
            show.signif=FALSE,
            main="Presence/Absence Blue PF (Every3)")

set.seed(7)
brt_PA_blue.every3.signif <- gbm.step(data = rvc_pf3, 
                                      gbm.x = c(13,16,17,18,19,21,24,25,33,34,37,39),
                                      gbm.y = 9,
                                      family = "bernoulli", 
                                      tree.complexity = 5,
                                      learning.rate = 0.01, 
                                      bag.fraction = 0.75)
Perf <- ggPerformance(blue.PA=brt_PA_blue.every3.signif) 
round(Perf,2)
ggInfluence(brt_PA_blue.every3.signif,col.bar = "#1380BE", 
            show.signif=FALSE,
            main="Presence/Absence Blue PF (Every3)")

brt_PA_blue.every3.signif.prerun<- plot.gbm.4list(brt_PA_blue.every3.signif)
brt_PA_blue.every3.signif.boot <- gbm.bootstrap.functions(brt_PA_blue.every3.signif, list.predictors=brt_PA_blue.every3.signif.prerun, n.reps=100)

# Moran's I for PA blue model every 3rd site #
distancesPA <- as.matrix(dist(cbind(rvc_pf3$Latitude, rvc_pf3$Longitude)))
distancesPA[1:5, 1:5]
distancesPA.inverse <- 1/distancesPA
diag(distancesPA.inverse) <- 0
distancesPA.inverse[is.infinite(distancesPA.inverse)] <- 0 
distancesPA.inverse[1:5, 1:5]
Moran.I(brt_PA_blue.every3.signif$residuals, distancesPA.inverse, na.rm=TRUE)





# BRT models: blue (Scarus coeruleus) Biomass ####
set.seed(12)
brt_bio_blue <- gbm.step(data = rvc_pf[rvc_pf$logp1_bluePF_biomass_kg>0,], 
                         gbm.x = pred.blue, 
                         gbm.y = 5, 
                         family = "gaussian", 
                         tree.complexity = 5,
                         learning.rate = 0.01, 
                         bag.fraction = 0.75) 

Perf <- ggPerformance(blue.biomass=brt_bio_blue) 
round(Perf,2)
ggInfluence(brt_bio_blue, col.bar = "#1380BE", 
            show.signif=FALSE,
            main="Biomass Blue PF")

# PA blue signif
set.seed(12)
brt_bio_blue.signif <- gbm.step(data = rvc_pf[rvc_pf$logp1_bluePF_biomass_kg>0,], 
                                gbm.x = c(16,20,22,25,37,39), 
                                gbm.y = 5, 
                                family = "gaussian", 
                                tree.complexity = 5,
                                learning.rate = 0.001, 
                                bag.fraction = 0.75)
Perf <- ggPerformance(blue.biomass=brt_bio_blue.signif) 
round(Perf,2)
ggInfluence(brt_bio_blue.signif,col.bar = "#1380BE", 
            show.signif=FALSE,
            main="Biomass Blue PF")

brt_bio_blue.signif.prerun<- plot.gbm.4list(brt_bio_blue.signif)
brt_bio_blue.signif.boot <- gbm.bootstrap.functions(brt_bio_blue.signif, list.predictors=brt_bio_blue.signif.prerun, n.reps=100)

# Moran's I for Biomass blue model
distances <- as.matrix(dist(cbind(rvc_pf[rvc_pf$logp1_bluePF_biomass_kg>0,]$Latitude, 
                                  rvc_pf[rvc_pf$logp1_bluePF_biomass_kg>0,]$Longitude)))
distances[1:5, 1:5]
distances.inverse <- 1/distances
diag(distances.inverse) <- 0
distances.inverse[1:5, 1:5]
Moran.I(brt_bio_blue.signif$residuals, distances.inverse, na.rm=TRUE)




# BRT models: rainbow (Scarus guacamaia) PA ####

set.seed(4)
brt_PA_rainbow <- gbm.step(data = rvc_pf, 
                           gbm.x = pred.rainbowPA, 
                           gbm.y = 10, 
                           family = "bernoulli", 
                           tree.complexity = 5,
                           learning.rate = 0.01, 
                           bag.fraction = 0.75) 

Perf <- ggPerformance(rainbow.PA=brt_PA_rainbow) 
round(Perf,2)
ggInfluence(brt_PA_rainbow, col.bar = "cadetblue3", 
            show.signif=FALSE,
            main="Presence/Absence Rainbow PF")

# PA rainbow signif
set.seed(4)
brt_PA_rainbow.signif <- gbm.step(data = rvc_pf, 
                                  gbm.x = c(16,17,18,19,20,21,23,24,33,40), 
                                  gbm.y = 10, 
                                  family = "bernoulli", 
                                  tree.complexity = 5,
                                  learning.rate = 0.01, 
                                  bag.fraction = 0.75)
Perf <- ggPerformance(rainbow.PA=brt_PA_rainbow.signif) 
round(Perf,2)
ggInfluence(brt_PA_rainbow.signif,col.bar = "cadetblue3", 
            show.signif=FALSE,
            main="Presence/Absence Rainbow PF")

brt_PA_rainbow.signif.prerun<- plot.gbm.4list(brt_PA_rainbow.signif)
brt_PA_rainbow.signif.boot <- gbm.bootstrap.functions(brt_PA_rainbow.signif, list.predictors=brt_PA_rainbow.signif.prerun, n.reps=100)

# Moran's I for Pres/Abs rainbow model
distancesPA <- as.matrix(dist(cbind(rvc_pf$Latitude, rvc_pf$Longitude)))
distancesPA[1:5, 1:5]
distancesPA.inverse <- 1/distancesPA
diag(distancesPA.inverse) <- 0
distancesPA.inverse[is.infinite(distancesPA.inverse)] <- 0 
distancesPA.inverse[1:5, 1:5]
Moran.I(brt_PA_rainbow.signif$residuals, distancesPA.inverse, na.rm=TRUE)




# BRT models: rainbow (Scarus guacamaia) Biomass ####
set.seed(4)
brt_bio_rainbow <- gbm.step(data = rvc_pf[rvc_pf$logp1_rainbowPF_biomass_kg>0,], 
                            gbm.x = pred.rainbow, 
                            gbm.y = 6, 
                            family = "gaussian", 
                            tree.complexity = 5,
                            learning.rate = 0.001, 
                            bag.fraction = 0.75) 

Perf <- ggPerformance(rainbow.biomass=brt_bio_rainbow) 
round(Perf,2)
ggInfluence(brt_bio_rainbow, col.bar = "cadetblue3", 
            show.signif=FALSE,
            main="Biomass Rainbow PF")

# PA rainbow signif
set.seed(4)
brt_bio_rainbow.signif <- gbm.step(data = rvc_pf[rvc_pf$logp1_rainbowPF_biomass_kg>0,], 
                                   gbm.x = c(12,13,17,19,20,21,22,23,24,32,34,35,37,40),
                                   gbm.y = 6, 
                                   family = "gaussian", 
                                   tree.complexity = 5,
                                   learning.rate = 0.001, 
                                   bag.fraction = 0.75)
Perf <- ggPerformance(rainbow.biomass=brt_bio_rainbow.signif) 
round(Perf,2)
ggInfluence(brt_bio_rainbow.signif,col.bar = "cadetblue3", 
            show.signif=FALSE,
            main="Biomass Rainbow PF")
gbm.plot(brt_bio_rainbow.signif, n.plots=10, write.title = FALSE)

brt_bio_rainbow.signif.prerun<- plot.gbm.4list(brt_bio_rainbow.signif)
brt_bio_rainbow.signif.boot <- gbm.bootstrap.functions(brt_bio_rainbow.signif, list.predictors=brt_bio_rainbow.signif.prerun, n.reps=100)

# Moran's I for Biomass rainbow model
distances <- as.matrix(dist(cbind(rvc_pf[rvc_pf$logp1_rainbowPF_biomass_kg>0,]$Latitude, 
                                  rvc_pf[rvc_pf$logp1_rainbowPF_biomass_kg>0,]$Longitude)))
distances[1:5, 1:5]
distances.inverse <- 1/distances
diag(distances.inverse) <- 0
distances.inverse[1:5, 1:5]
Moran.I(brt_bio_rainbow.signif$residuals, distances.inverse, na.rm=TRUE)







# BRT models: stoplight (Sparisoma viride) PA ####

set.seed(12)
brt_PA_stoplight <- gbm.step(data = rvc_pf, 
                             gbm.x = pred.stoplightPA, 
                             gbm.y = 11, 
                             family = "bernoulli", 
                             tree.complexity = 5,
                             learning.rate = 0.01, 
                             bag.fraction = 0.75) 

Perf <- ggPerformance(stoplight.PA=brt_PA_stoplight) 
round(Perf,2)
ggInfluence(brt_PA_stoplight, col.bar = "#54D5B1", 
            show.signif=FALSE,
            main="Presence/Absence Stoplight PF")

# PA stoplight signif
set.seed(12)
brt_PA_stoplight.signif <- gbm.step(data = rvc_pf, 
                                    gbm.x = c(13,16,17,18,19,20,21,23,24,25,28,32,33,34,37,41), 
                                    gbm.y = 11, 
                                    family = "bernoulli", 
                                    tree.complexity = 5,
                                    learning.rate = 0.01, 
                                    bag.fraction = 0.75)
Perf <- ggPerformance(stoplight.PA=brt_PA_stoplight.signif) 
round(Perf,2)
ggInfluence(brt_PA_stoplight.signif,col.bar = "#54D5B1", 
            show.signif=FALSE,
            main="Presence/Absence Stoplight PF")

brt_PA_stoplight.signif.prerun<- plot.gbm.4list(brt_PA_stoplight.signif)
brt_PA_stoplight.signif.boot <- gbm.bootstrap.functions(brt_PA_stoplight.signif, list.predictors=brt_PA_stoplight.signif.prerun, n.reps=100)

# Moran's I for Pres/Abs stoplight model
distancesPA <- as.matrix(dist(cbind(rvc_pf$Latitude, rvc_pf$Longitude)))
distancesPA[1:5, 1:5]
distancesPA.inverse <- 1/distancesPA
diag(distancesPA.inverse) <- 0
distancesPA.inverse[is.infinite(distancesPA.inverse)] <- 0 
distancesPA.inverse[1:5, 1:5]
Moran.I(brt_PA_stoplight.signif$residuals, distancesPA.inverse, na.rm=TRUE)

# BRT models: stoplight (Sparisoma viride) Biomass ####

set.seed(8)
brt_bio_stoplight <- gbm.step(data = rvc_pf[rvc_pf$logp1_stoplightPF_biomass_kg>0,], 
                              gbm.x = pred.stoplight, 
                              gbm.y = 7, 
                              family = "gaussian", 
                              tree.complexity = 5,
                              learning.rate = 0.001, 
                              bag.fraction = 0.75) 

Perf <- ggPerformance(stoplight.biomass=brt_bio_stoplight) 
round(Perf,2)
ggInfluence(brt_bio_stoplight, col.bar = "#54D5B1", 
            show.signif=FALSE,
            main="Biomass Stoplight PF")

# PA stoplight signif
set.seed(8)
brt_bio_stoplight.signif <- gbm.step(data = rvc_pf[rvc_pf$logp1_stoplightPF_biomass_kg>0,], 
                                     gbm.x = c(12,13,16,17,18,19,20,21,22,23,24,25,28,33,34,36,37,41), 
                                     gbm.y = 7, 
                                     family = "gaussian", 
                                     tree.complexity = 5,
                                     learning.rate = 0.001, 
                                     bag.fraction = 0.75)
Perf <- ggPerformance(stoplight.biomass=brt_bio_stoplight.signif) 
round(Perf,2)
ggInfluence(brt_bio_stoplight.signif,col.bar = "#54D5B1", 
            show.signif=FALSE,
            main="Biomass Stoplight PF")

brt_bio_stoplight.signif.prerun<- plot.gbm.4list(brt_bio_stoplight.signif)
brt_bio_stoplight.signif.boot <- gbm.bootstrap.functions(brt_bio_stoplight.signif, list.predictors=brt_bio_stoplight.signif.prerun, n.reps=100)

# Moran's I for Biomass stoplight model
distances <- as.matrix(dist(cbind(rvc_pf[rvc_pf$logp1_stoplightPF_biomass_kg>0,]$Latitude, 
                                  rvc_pf[rvc_pf$logp1_stoplightPF_biomass_kg>0,]$Longitude)))
distances[1:5, 1:5]
distances.inverse <- 1/distances
diag(distances.inverse) <- 0
distances.inverse[1:5, 1:5]
distances.inverse[is.infinite(distances.inverse)] <- 0 
Moran.I(brt_bio_stoplight.signif$residuals, distances.inverse, na.rm=TRUE)






# Comparison table plots ####

# Change a few things about the ggMultiInfluence function:
# make the background on the plot white rather than grey50 (in the scale_fill_gradient na.value= call)
# specify the order of variables and change the names

ggMultiInfluence1<-function (...,col.gradient=c("white","lightblue4"),round=1,
                             col.text="grey10",col.grid=NULL,size.grid=0.3,legend.pos="bottom",
                             legend.dir="horizontal", scale.gradient=c(0,100)){
  
  BRT<-list(...)
  nBRT<-length(BRT)
  dfContr<-list()
  
  # Extract the relative influence of each predictor for each model
  for (i in 1:nBRT){
    dfContr[[i]]<-data.frame(BRT[[i]]$contributions)
    if(is.null(names(BRT))){
      colnames(dfContr[[i]])<-c("Predictor",paste("Model", i, sep = " "))}
    else{
      colnames(dfContr[[i]])<-c("Predictor",names(BRT)[i])}
  }
  
  # Merge them into a single dataframe
  Influence<-join_all(dfContr, by = 'Predictor')
  
  # Transform the data frame and round the values
  I.melted <- melt(Influence)
  I.melted$value<-round(I.melted$value,round)
  
  I.sorted <- I.melted %>%
    arrange(Predictor) %>%
    mutate(Predictor = factor(Predictor, levels=c(
      "random",
      "Year",
      "Month",
      "RecEng",
      "Marine_reserve",
      "Grav_tot_closest",
      "lspop_20km",
      "wave_exposure",
      "ALL_TURF_200",
      "SST",
      "Reef_complexity",
      "pf_preds_biomass_kg",
      "pf_no_midnight_biomass_kg",
      "nursery_seagrass",
      "nursery_mangroves",
      "NPP",
      "ALL_MACROALGAE_200",
      "Habitat_type_classLV2",
      "Depth",
      "Deepwater",
      "Coral_cover",
      "Coral_area_UFRTM_20km",
      "Coral_area_UFRTM_200km",
      "connectivity",
      "Artificial_reefs_1km")))
  
  labels <- c("Random",
              "Year",
              "Month",
              "Rec fishing engagement",
              "Marine reserve",
              "Market gravity",
              "Human population 20km",
              "Wave exposure",
              "Turf cover",
              "SST",
              "Reef complexity",
              "Parrotfish predators",
              "Other parrotfish",
              "Nursery availability (seagrass)",
              "Nursery availability (mangroves)",
              "NPP",
              "Macroalgae cover",
              "Habitat type",
              "Depth",
              "Deepwater",
              "Coral cover",
              "Coral area 20km",
              "Coral area 200km",
              "Connectivity",
              "Artificial reefs")
  
  # Plot the heatmap
  p<-ggplot(I.sorted, aes(x = variable, y = Predictor))+
    geom_tile(aes(fill=value),size=size.grid)+
    geom_text(aes(label=value,size=value), color=col.text,show.legend = FALSE)+
    scale_x_discrete(position = "top")+
    scale_y_discrete(labels=labels)+
    scale_fill_gradient(low=col.gradient[1],high=col.gradient[2],name="Relative influence (%)",
                        limits=scale.gradient, na.value="white")+
    labs(x= "",y = "")+
    theme_minimal()+
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(face="bold"),
          legend.position=legend.pos,legend.direction=legend.dir)
  
  if(!is.null(col.grid)){
    p<-p+geom_tile(aes(fill=value),colour=col.grid,size=size.grid)+
      geom_text(aes(label=value,size=value), color=col.text,show.legend = FALSE)}
  
  print(p)
  print(Influence,digits = 3)
}


# Plot relative influence of predictors for three single species pres/abs models
# Dummary is needed so that all variables show up on the left axis. The dummy column is then edited out
# in Powerpoint or Illustrator. 

ggMultiInfluence1(Dummy=brt_PA_midnight,
                  Midnight=brt_PA_midnight.signif,
                  Blue=brt_PA_blue.signif,
                  Rainbow=brt_PA_rainbow.signif,
                  Stoplight=brt_PA_stoplight.signif,
                  col.gradient=c("white", "#011a42"),
                  col.text="black")

# Plot relative influence of predictors for three single species biomass models
ggMultiInfluence1(Dummy=brt_bio_midnight,
                  Midnight=brt_bio_midnight.signif,
                  Blue=brt_bio_blue.signif,
                  Rainbow=brt_bio_rainbow.signif,
                  Stoplight=brt_bio_stoplight.signif,
                  col.gradient=c("white", "#011a42"),
                  col.text="black")

# Plot relative influence of predictors for the original and the 'every 3rd site by latitude' dataset
ggMultiInfluence1(Dummy=brt_PA_midnight,
                  Midnight=brt_PA_midnight.signif, 
                  Midnight_subsample=brt_PA_midnight.every3.signif,
                  Blue=brt_PA_blue.signif, 
                  Blue_subsample=brt_PA_blue.every3.signif,
                  col.gradient=c("white", "#011a42"),
                  col.text="black")









# Partial dependency plots ####

# ggPD_boot_ksrz: function for plotting the bootstrap

# amended from J. Jouffrey's original ggPD_boot function: 
#https://github.com/JBjouffray/ggBRT/blob/master/R/ggPD_boot.R

# Changes made by:
# Katie Sievers, James Cook University (code provided by Eva McClure),
# Rachel Zuercher, Florida International University (rachel.zuercher@gmail.com)

# Katie Sievers changes to ggPD_boot():
# KS change to add *rug.col* (line 849) so that you can change the color of the rugs
# (default is set to black)
# KS removed ,... after TRUE on line 895

# Rachel Zuercher changes to ggPD_boot():
# RZ added 'na.rm=T' within the mean() function on lines 905 and 908 (to make the temp dataframe) 
# for cases where the booted.preds object has NAs
# RZ edited ggplot() code lines 952-964 and added 970-972 to plot rugs and confidence intervals on variables 
# that are factors (where no predictor has been defined in the ggPD_boot() function)
# RZ added lines 1019-1021 to specify plotting instructions if rug=T for is.null(predictor)=TRUE

ggPD_boot_ksrz<- function(gbm.object, predictor = NULL,n.plots = length(pred.names),list.4.preds=NULL,booted.preds=NULL, nrow=NULL,ncol=NULL,
                          col.line="darkorange",cex.line=0.5, type.ci="lines",col.ci= "grey80",cex.ci=0.3,lty.ci=2, alpha.ci=0.1,smooth = FALSE,
                          col.smooth="blue",cex.smooth=0.3,span=0.3,rug = FALSE, rug.col = 'grey80', rug.pos="t",common.scale = TRUE,cis=c(0.025, 0.975),
                          y.label = "Fitted function",x.label=NULL,...){
  
  gbm.call <- gbm.object$gbm.call
  pred.names <- gbm.call$predictor.names
  
  ggPD_boot.plots<-function (gbm.object) {
    if (!requireNamespace("gbm")) {
      stop("you need to install the gbm package to run this function")
    }
    
    if (is.null(booted.preds)) {
      stop("you need to set booted.preds as the array from the bootstrap run
           (eg testboot$function.preds using testboot<-gbm.bootstrap.functions())")
    }
    
    if (is.null(list.4.preds)) {
      stop("you need to set list.4.preds as the result of plot.gbm.4list()")
    }
    requireNamespace("splines")
    gbm.x <- gbm.call$gbm.x
    response.name <- gbm.call$response.name
    nt<-gbm.object$n.trees
    data <- gbm.call$dataframe  #data set using for analysis (fish data enviro)
    
    max.vars <- length(gbm.object$contributions$var)
    if (n.plots > max.vars) {
      n.plots <- max.vars
      warning("reducing no of plotted predictors to maximum available (",
              max.vars, ")")
    }
    predictors <- list(rep(NA, n.plots))
    responses <- list(rep(NA, n.plots))
    responses.lower <- list(rep(NA, n.plots))
    responses.upper <- list(rep(NA, n.plots))
    
    for (j in c(1:max.vars)) {
      k <- match(gbm.object$contributions$var[j], pred.names)
      
      if (is.null(x.label)) {
        var.name <- gbm.call$predictor.names[k]
      }
      else {
        var.name <- x.label
      }
      pred.data <- data[, gbm.call$gbm.x[k]]
      response.matrix <- gbm::plot.gbm(gbm.object, i.var = k, n.trees = nt, return.grid = TRUE) #,... took out at end after true
      predictors[[j]] <- response.matrix[, 1]
      if (is.factor(data[, gbm.call$gbm.x[k]])) {
        predictors[[j]] <- factor(predictors[[j]], levels = levels(data[,gbm.call$gbm.x[k]]))
      }
      
      responses[[j]] <- response.matrix[, 2] - mean(response.matrix[, 2])
      
      num.values <- nrow(response.matrix)
      
      temp <-apply(booted.preds[,k,] - mean(booted.preds[,k,], na.rm=T), 1, function(x){quantile(x, cis[1],na.rm=T)})
      responses.lower[[j]] <- temp[1:num.values]
      
      temp <-apply(booted.preds[,k,] - mean(booted.preds[,k,], na.rm=T), 1, function(x){quantile(x, cis[2],na.rm=T)})
      responses.upper[[j]] <- temp[1:num.values]
      
      if(j == 1) {
        ymin = min(responses.lower[[j]])
        ymax = max(responses.upper[[j]])
        dat<-data.frame(pred.data)
      }
      else {
        ymin = min(ymin,min(responses.lower[[j]]))
        ymax = max(ymax,max(responses.upper[[j]]))
        dat<-data.frame(dat,pred.data)
      }
    }
    
    if (is.null(predictor)){
      
      fittedFunc<-list()
      fittedFunc.lower<-list()
      fittedFunc.upper<-list()
      fittedVal<-list()
      ribbon<-list()
      ggPD<-list()
      
      for (i in 1:n.plots) {
        k <- match(gbm.object$contributions$var[i], pred.names)
        var.name <- gbm.call$predictor.names[k]
        
        fittedFunc[[i]]<-data.frame(predictors[i],responses[i])
        colnames(fittedFunc[[i]])<-c("x","y")
        
        fittedFunc.lower[[i]]<-data.frame(predictors[i],responses.lower[i])
        colnames(fittedFunc.lower[[i]])<-c("x","y")
        
        fittedFunc.upper[[i]]<-data.frame(predictors[i],responses.upper[i])
        colnames(fittedFunc.upper[[i]])<-c("x","y")
        
        fittedVal[[i]]<-data.frame(gbm.object$fitted,dat[i])
        colnames(fittedVal[[i]])<-c("y","x")
        
        ribbon[[i]]<-data.frame("x"=fittedFunc.lower[[i]]$x,"ylow"=fittedFunc.lower[[i]]$y,"yup"=fittedFunc.upper[[i]]$y)
        
        if (is.factor(fittedFunc[[i]]$x)) {
          
          ggPD[[i]]<- ggplot()+
            geom_point(data=fittedFunc[[i]], aes(x=x,y=y), col=col.line, size=2)+
            geom_errorbar(data=ribbon[[i]], aes(x=x, ymin=ylow, ymax=yup), width=.2,
                          col=col.line, alpha=0.4)+
            ylab(y.label)+
            xlab(paste(var.name, "  (", round(gbm.object$contributions[i,2], 1), "%)", sep = ""))+
            theme_bw()+
            theme(panel.grid.minor = element_line(linetype = "blank"),
                  panel.grid.major = element_line(linetype = "blank"),
                  axis.text.x  = element_text(size=6),
                  axis.title.x  = element_text(size=10),
                  axis.line.y = element_line(size=0.1),
                  axis.line.x=element_line(size=0.1))
          
          if (common.scale==T){
            ggPD[[i]]<-ggPD[[i]]+ylim(c(ymin,ymax))
          }
          
          if (rug==T){
            ggPD[[i]]<-ggPD[[i]]+geom_rug(data=fittedVal[[i]],aes(x=x,y=y),sides=rug.pos,position="jitter",color=rug.col)
          }
        }
        
        if (type.ci=="lines"){
          ggPD[[i]]<- ggplot(fittedFunc[[i]], aes(x=x,y=y))+
            geom_line(color=col.line,size=cex.line)+
            geom_line(data=fittedFunc.lower[[i]],aes(x=x,y=y),size=cex.ci,color=col.ci,linetype=lty.ci)+
            geom_line(data=fittedFunc.upper[[i]],aes(x=x,y=y),size=cex.ci,color=col.ci,linetype=lty.ci)+
            ylab(y.label)+
            xlab(paste(var.name, "  (", round(gbm.object$contributions[i,2], 1), "%)", sep = ""))+
            theme_bw()+
            theme(panel.grid.minor = element_line(linetype = "blank"),
                  panel.grid.major = element_line(linetype = "blank"),
                  axis.title.x  = element_text(size=10),
                  axis.line.y = element_line(size=0.1),
                  axis.line.x=element_line(size=0.1))
          
          if (smooth==T){
            ggPD[[i]]<-ggPD[[i]]+geom_smooth(span=span,size=0.3,color=col.smooth,se=F,linetype=2)
          }
          
          if (rug==T){
            ggPD[[i]]<-ggPD[[i]]+geom_rug(data=fittedVal[[i]],aes(x=x,y=y),sides=rug.pos,position="jitter",color=rug.col)
          }
          
          if (common.scale==T){
            ggPD[[i]]<-ggPD[[i]]+ylim(c(ymin,ymax))
          }
        }
        
        if (type.ci=="ribbon" & is.numeric(fittedFunc[[i]]$x)) {
          ggPD[[i]]<- ggplot()+
            geom_ribbon(data=ribbon[[i]],aes(x=x,ymin=ylow,ymax=yup),fill=col.ci,alpha=alpha.ci)+
            geom_line(data=fittedFunc[[i]], aes(x=x,y=y),color=col.line,size=cex.line)+
            ylab(y.label)+
            xlab(paste(var.name, "  (", round(gbm.object$contributions[i,2], 1), "%)", sep = ""))+
            theme_bw()+
            theme(panel.grid.minor = element_line(linetype = "blank"),
                  panel.grid.major = element_line(linetype = "blank"),
                  axis.title.x  = element_text(size=10),
                  axis.line.y = element_line(size=0.1),
                  axis.line.x=element_line(size=0.1))
          
          if (smooth==T){
            ggPD[[i]]<-ggPD[[i]]+geom_smooth(data=fittedFunc[[i]],aes(x=x,y=y),span=span,size=0.3,color=col.smooth,se=F,linetype=2)
          }
          
          if (rug==T){
            ggPD[[i]]<-ggPD[[i]]+geom_rug(data=fittedVal[[i]],aes(x=x,y=y),sides=rug.pos,position="jitter",color=rug.col)
          }
          
          if (common.scale==T){
            ggPD[[i]]<-ggPD[[i]]+ylim(c(ymin,ymax))
          }
        }
      }
      list(ggPD=ggPD)
    }
    
    else{
      
      if (is.character(predictor)){
        predictor<-match(predictor,gbm.object$contributions$var)}
      
      k <- match(gbm.object$contributions$var[predictor], pred.names)
      var.name <- gbm.call$predictor.names[k]
      
      fittedFunc<-data.frame(predictors[predictor],responses[predictor])
      colnames(fittedFunc)<-c("x","y")
      
      fittedFunc.lower<-data.frame(predictors[predictor],responses.lower[predictor])
      colnames(fittedFunc.lower)<-c("x","y")
      
      fittedFunc.upper<-data.frame(predictors[predictor],responses.upper[predictor])
      colnames(fittedFunc.upper)<-c("x","y")
      
      ribbon<-data.frame("x"=fittedFunc.lower$x,"ylow"=fittedFunc.lower$y,"yup"=fittedFunc.upper$y)
      
      fittedVal<-data.frame(gbm.object$fitted,dat[predictor])
      colnames(fittedVal)<-c("y","x")
      
      if (is.factor(fittedFunc$x)) {
        ggPD<- ggplot(fittedFunc, aes(x=x,y=y))+
          geom_boxplot(color=col.line,size=cex.line)+
          geom_boxplot(data=fittedFunc.lower, aes(x=x,y=y),color=col.ci)+
          geom_boxplot(data=fittedFunc.upper, aes(x=x,y=y),color=col.ci)+
          ylab(y.label)+
          xlab(paste(var.name, "  (", round(gbm.object$contributions[predictor,2], 1), "%)", sep = ""))+
          theme_bw()+
          theme(panel.grid.minor = element_line(linetype = "blank"),
                panel.grid.major = element_line(linetype = "blank"),
                axis.text.x  = element_text(size=6),
                axis.title.x  = element_text(size=10),
                axis.line.y = element_line(size=0.1),
                axis.line.x=element_line(size=0.1))
        
        if (common.scale==T){
          ggPD<-ggPD+ylim(c(ymin,ymax))}
        
        if (rug==T){
          ggPD<-ggPD+geom_rug(data=fittedVal,aes(x=x,y=y),sides=rug.pos,position="jitter",color=rug.col)
        }
      }
      
      if (type.ci=="lines"){
        ggPD<- ggplot(fittedFunc, aes(x=x,y=y))+
          geom_line(color=col.line,size=cex.line)+
          geom_line(data=fittedFunc.lower,aes(x=x,y=y),size=cex.ci,color=col.ci,linetype=lty.ci)+
          geom_line(data=fittedFunc.upper,aes(x=x,y=y),size=cex.ci,color=col.ci,linetype=lty.ci)+
          ylab(y.label)+
          xlab(paste(var.name, "  (", round(gbm.object$contributions[predictor,2], 1), "%)", sep = ""))+
          theme_bw()+
          theme(panel.grid.minor = element_line(linetype = "blank"),
                panel.grid.major = element_line(linetype = "blank"),
                axis.title.x  = element_text(size=10),
                axis.line.y = element_line(size=0.1),
                axis.line.x=element_line(size=0.1))
        
        if (smooth==T){
          ggPD<-ggPD+geom_smooth(span=span,size=0.3,color=col.smooth,se=F,linetype=2)
        }
        
        if (rug==T){
          ggPD<-ggPD+geom_rug(data=fittedVal,aes(x=x,y=y),sides=rug.pos,position="jitter",color=rug.col)
        }
        
        if (common.scale==T){
          ggPD<-ggPD+ylim(c(ymin,ymax))
        }
      }
      
      if (type.ci=="ribbon"){
        ggPD<- ggplot()+
          geom_ribbon(data=ribbon,aes(x=x,ymin=ylow,ymax=yup),fill=col.ci,alpha=alpha.ci)+
          geom_line(data=fittedFunc,aes(x=x,y=y),color=col.line,size=cex.line)+
          ylab(y.label)+
          xlab(paste(var.name, "  (", round(gbm.object$contributions[predictor,2], 1), "%)", sep = ""))+
          theme_bw()+
          theme(panel.grid.minor = element_line(linetype = "blank"),
                panel.grid.major = element_line(linetype = "blank"),
                axis.title.x  = element_text(size=10),
                axis.line.y = element_line(size=0.1),
                axis.line.x=element_line(size=0.1))
        
        if (smooth==T){
          ggPD<-ggPD+geom_smooth(data=fittedFunc,aes(x=x,y=y),span=span,size=0.3,color=col.smooth,se=F,linetype=2)
        }
        
        if (rug==T){
          ggPD<-ggPD+geom_rug(data=fittedVal,aes(x=x,y=y),sides=rug.pos,position="jitter",color=rug.col)
        }
        
        if (common.scale==T){
          ggPD<-ggPD+ylim(c(ymin,ymax))
        }
      }
      list(ggPD=ggPD)
    }
    }
  
  plot<-ggPD_boot.plots(gbm.object)
  plot
  
  if(is.null(predictor)){
    do.call(grid.arrange,c(plot$ggPD,list(nrow=nrow,ncol=ncol)))}
  else grid.draw(plot$ggPD)
}


# Midnight PA
plot.midnightPA <- ggPD_boot_ksrz(brt_PA_midnight.signif, list.4.preds=brt_PA_midnight.signif.prerun, 
                                  booted.preds=brt_PA_midnight.signif.boot$function.preds, 
                                  n.plots = 9, nrow=2, ncol=5, 
                                  col.line="midnightblue", cex.line=1,  
                                  cis =c(0.05, 0.95), cex.ci =1, col.ci = "midnightblue", 
                                  alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                                  rug = TRUE, rug.col = 'darkgrey',
                                  y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("midnightPA.pdf", plot.midnightPA, units="mm", width=290, height=100, device=cairo_pdf)

# Blue PA
plot.bluePA <- ggPD_boot_ksrz(brt_PA_blue.signif, list.4.preds=brt_PA_blue.signif.prerun, 
                              booted.preds=brt_PA_blue.signif.boot$function.preds, 
                              n.plots = 7, nrow=2, ncol=5, 
                              col.line="#1380BE", cex.line=1,  
                              cis =c(0.05, 0.95), cex.ci =1, col.ci = "#1380BE", 
                              alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                              rug = TRUE, rug.col = 'darkgrey',
                              y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("bluePA.pdf", plot.bluePA, units="mm", width=290, height=100, device=cairo_pdf)

# Rainbow PA
plot.rainbowPA <- ggPD_boot_ksrz(brt_PA_rainbow.signif, list.4.preds=brt_PA_rainbow.signif.prerun, 
                                 booted.preds=brt_PA_rainbow.signif.boot$function.preds, 
                                 n.plots = 7, nrow=2, ncol=5, 
                                 col.line="cadetblue3", cex.line=1,  
                                 cis =c(0.05, 0.95), cex.ci =1, col.ci = "cadetblue3", 
                                 alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                                 rug = TRUE, rug.col = 'darkgrey',
                                 y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("rainbowPA.pdf", plot.rainbowPA, units="mm", width=290, height=100, device=cairo_pdf)

# Stoplight PA
plot.stoplightPA <- ggPD_boot_ksrz(brt_PA_stoplight.signif, list.4.preds=brt_PA_stoplight.signif.prerun, 
                                   booted.preds=brt_PA_stoplight.signif.boot$function.preds, 
                                   n.plots = 4, nrow=2, ncol=5, 
                                   col.line="#54D5B1", cex.line=1,  
                                   cis =c(0.05, 0.95), cex.ci =1, col.ci = "#54D5B1", 
                                   alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                                   rug = TRUE, rug.col = 'darkgrey',
                                   y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("stoplightPA.pdf", plot.stoplightPA, units="mm", width=290, height=100, device=cairo_pdf)






# Midnight Biomass
plot.midnightBio <- ggPD_boot_ksrz(brt_bio_midnight.signif, list.4.preds=brt_bio_midnight.signif.prerun, 
                                   booted.preds=brt_bio_midnight.signif.boot$function.preds, 
                                   n.plots = 3, nrow=2, ncol=5, 
                                   col.line="midnightblue", cex.line=1,  
                                   cis =c(0.05, 0.95), cex.ci =1, col.ci = "midnightblue", 
                                   alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                                   rug = TRUE, rug.col = 'darkgrey',
                                   y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("midnightBio.pdf", plot.midnightBio, units="mm", width=290, height=100, device=cairo_pdf)

# Blue Biomass
plot.blueBio <- ggPD_boot_ksrz(brt_bio_blue.signif, list.4.preds=brt_bio_blue.signif.prerun, 
                               booted.preds=brt_bio_blue.signif.boot$function.preds, 
                               n.plots = 6, nrow=2, ncol=5, 
                               col.line="#1380BE", cex.line=1,  
                               cis =c(0.05, 0.95), cex.ci =1, col.ci = "#1380BE", 
                               alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                               rug = TRUE, rug.col = 'darkgrey',
                               y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("blueBio.pdf", plot.blueBio, units="mm", width=290, height=100, device=cairo_pdf)

# Rainbow Biomass
plot.rainbowBio <- ggPD_boot_ksrz(brt_bio_rainbow.signif, list.4.preds=brt_bio_rainbow.signif.prerun, 
                                  booted.preds=brt_bio_rainbow.signif.boot$function.preds, 
                                  n.plots = 10, nrow=2, ncol=5, 
                                  col.line="cadetblue3", cex.line=1,  
                                  cis =c(0.05, 0.95), cex.ci =1, col.ci = "cadetblue3", 
                                  alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                                  rug = TRUE, rug.col = 'darkgrey',
                                  y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("rainbowBio.pdf", plot.rainbowBio, units="mm", width=290, height=100, device=cairo_pdf)

# Stoplight Biomass
plot.stoplightBio <- ggPD_boot_ksrz(brt_bio_stoplight.signif, list.4.preds=brt_bio_stoplight.signif.prerun, 
                                    booted.preds=brt_bio_stoplight.signif.boot$function.preds, 
                                    n.plots = 6, nrow=2, ncol=5, 
                                    col.line="#54D5B1", cex.line=1,  
                                    cis =c(0.05, 0.95), cex.ci =1, col.ci = "#54D5B1", 
                                    alpha.dot = 0.2, type.ci = 'ribbon', alpha.ci = 0.1,
                                    rug = TRUE, rug.col = 'darkgrey',
                                    y.label = "", common.scale = TRUE, smooth = FALSE) 
ggsave("stoplightBio.pdf", plot.stoplightBio, units="mm", width=290, height=100, device=cairo_pdf)




