#########################################
### Final variable selection and model ##
#########################################

.libPaths("D:/R libraries")
library(lme4)
library(MuMIn)
library(visreg)
library(effects)
library(dplyr)
library(jtools)
library(ggplot2)
library(ggeffects)
require(GGally)
require(reshape2)
require(compiler)
require(parallel)
require(boot)
require(lattice)
require(stats)
library(sjmisc)
library(plyr)
library(mgcv)

##### data categorization and scaled #####

pkcovar <- read.csv("pk_original_EJ.csv")

# macaques: mf vs other

pkcovar$host_group <- revalue(pkcovar$host_species, c("Macaca nemestrina"="Other", "Macaca arctoides"="Other", "Semnopithecus obscurus"="Other", "Macaca leonina"="Other",
                                                      "Macaca fascicularis"="Macaca fascicularis"))
pkcovar$host_group <- factor(pkcovar$host_group)
pkcovar$host_group <- factor(pkcovar$host_group, c("Other", "Macaca fascicularis"))

# scale elevation and population

pkcovar$scale5km <- scale(pkcovar$X5km_srtm_mean, center = TRUE, scale = TRUE)
pkcovar$scale10km <- scale(pkcovar$X10km_srtm_mean, center = TRUE, scale = TRUE)
pkcovar$scale20km <- scale(pkcovar$X20km_srtm_mean, center = TRUE, scale = TRUE)

pkcovar$scale_pop5km <- scale(pkcovar$X5km_pop12_mean, center = TRUE, scale = TRUE)
pkcovar$scale_pop10km <- scale(pkcovar$X10km_pop12_mean, center = TRUE, scale = TRUE)
pkcovar$scale_pop20km <- scale(pkcovar$X20km_pop12_mean, center = TRUE, scale = TRUE)

# forest 

pkcovar$forest5km <- scale(pkcovar$fc5km, center = TRUE, scale = TRUE)
pkcovar$forest10km <- scale(pkcovar$fc10km, center = TRUE, scale = TRUE)
pkcovar$forest20km <- scale(pkcovar$fc20km, center = TRUE, scale = TRUE)

pkcovar$fc_5km_discr.   <- cut(pkcovar$fc5km, 
                               breaks=c(-Inf, 0.2, 0.5, Inf), 
                               labels=c("low","med", "high"))
pkcovar$fc_5km_discr. <- factor(pkcovar$fc_5km_discr.)

pkcovar$fc_10km_discr.  <- cut(pkcovar$fc10km, 
                               breaks=c(-Inf, 0.2, 0.5, Inf), 
                               labels=c("low","med", "high"))
pkcovar$fc_10km_discr. <- factor(pkcovar$fc_10km_discr.)

pkcovar$fc_20km_discr.  <- cut(pkcovar$fc20km, 
                               breaks=c(-Inf, 0.2, 0.5, Inf), 
                               labels=c("low","med", "high"))
pkcovar$fc_20km_discr. <- factor(pkcovar$fc_20km_discr.)

# fragmentation

pkcovar$frag5km <- scale(pkcovar$para_5km, center = TRUE, scale = TRUE)
pkcovar$frag10km <- scale(pkcovar$para_10km, center = TRUE, scale = TRUE)
pkcovar$frag20km <- scale(pkcovar$para_20km, center = TRUE, scale = TRUE)

pkcovar$para5_q <- with(pkcovar, factor(
  findInterval(para5km, c(-Inf,
                          quantile(para5km, probs=c(0.25, .5, .75)), Inf)), 
  labels=c("Q1","Q2","Q3","Q4")
))

pkcovar$para10_q <- with(pkcovar, factor(
  findInterval(para10km, c(-Inf,
                           quantile(para10km, probs=c(0.25, .5, .75)), Inf)), 
  labels=c("Q1","Q2","Q3","Q4")
))

pkcovar$para20_q <- with(pkcovar, factor(
  findInterval(para20km, c(-Inf,
                           quantile(para20km, probs=c(0.25, .5, .75)), Inf)), 
  labels=c("Q1","Q2","Q3","Q4")
))

####################################
# define function for VIF stepwise #
####################################

vif_func<-function(in_frame,thresh,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
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
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
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

##### full model, scaled variables #####

f0 <- glmer(cbind(pk_all, denom_tested-pk_all) ~ host_group + scale(srtm5km_mean) + scale(srtm10km_mean) + scale(srtm20km_mean) +
              scale(pop5km_mean) + scale(pop10km_mean) + scale(pop20km_mean) +
              scale(para_5km) + scale(para_10km) +scale(para_20km) +
              I((scale(para_5km))^2) + I((scale(para_10km))^2) + I((scale(para_20km))^2) +
              scale(fc5km) + scale(fc10km) + scale(fc20km) +
              (1|source_id),
            data=pkcovar10,
            family="binomial")

# VIF < 6

variables_6<-vif_func(in_frame=pkcorr, thresh=6)
variables_6

############
### AIC ####
############

# Forwards selection strategy

mod.0 <- glmer(cbind(pk_all, denom_tested-pk_all) ~ 1 + (1|source_id),
               data=pkcovar10, family="binomial")

mod.f <- glmer(cbind(pk_all, denom_tested-pk_all) ~ 1 + 
                 fc20km + 
              #   para_5km+I(para_5km^2) + 
              #   host_group +  
              #   para_20km+I(para_20km^2) + 
              #   pop5km_mean + 
              #   fc5km + 
              #   srtm5km_mean + 
                 (1|source_id),
               data=pkcovar10, family="binomial")

summary(mod.f)
anova(mod.0, mod.f)

est <- coef(summary(mod.f))


###################
### final model ###
###################


mod.vif_6<- glmer(cbind(pk_all, denom_tested-pk_all) ~
                    srtm5km_mean+
                    pop5km_mean+
                    fc5km+
                    fc20km+
                    para_5km+I(para_5km^2)+
                    para_20km+I(para_20km^2)+
                    (1|source_id)+
                    host_group,
                  data=pkcovar10,
                  family="binomial")


summary(mod.vif_6)
est <- coef(summary(mod.vif_6))


######################
### visualization ####
######################

# plot linear variables 
m <- plot_model(mod.vif_6,
                terms = c("host_group [Macaca fascicularis]", "srtm5km_mean", "pop5km_mean", "fc5km", "fc20km"),
                title = "Adjusted OR and 95% CI for P. knowlesi linear risk factors at influential scales                                                    ",
                axis.lim = c(0.5, 9),
                axis.labels = c("M. fascicularis", "Forest [20km]", "Forest [5km]", "Population [5km]", "Elevation [5km]"), 
                wrap.title = 100, 
                dot.size = 4, 
                colors = c("seagreen", "navy"),
                vline.color = "red",
                width = .3) 
m

# plot deforestation at 20km with ggplot
library(predict3d)

pred_20km <- as.data.frame(ggpredict(mod.vif_6, term="para_20km [all]"))
head(pred_20km) # dataset with the prediciton for the plot

ggplot(pred_20km, aes(x = x, y = predicted)) + 
  geom_line(col= 'navy', size = 1) + 
  geom_ribbon(aes(ymin =  conf.low, ymax =  conf.high), alpha = 0.1, fill='navy') +
  theme_classic() +
  labs(x = "P:A ratio [20km]", y = "Prediction")



