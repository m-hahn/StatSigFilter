## ----setup,include=FALSE,cache=FALSE,echo=FALSE--------------------------
library(knitr)
library(coda)
library(plyr)
library(ggplot2)
library(xtable)
library(dplyr)
library(SIN)
library(papaja)

library(rstan)
library(parallel)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

library(loo)
library(lme4)

library(rstantools)

opts_chunk$set(fig.path='figures/fig-', 
               fig.align='center', 
               fig.show='hold',warning=FALSE,include=FALSE,message=FALSE)
options(replace.assign=TRUE,show.signif.stars=FALSE)
options(replace.assign=TRUE,width=75)
#opts_chunk$set(dev='postscript')

library(brms)

## ----loadfunctions, include=TRUE, echo=FALSE, warning=FALSE, error=TRUE, message=TRUE----
source("../R/createStanDat.R")
source("../R/createStanDatAcc.R")
source("../R/magnifytext.R")
source("../R/multiplot.R")
source("../R/plotresults.R")
source("../R/plotpredictions.R")
source("../R/stan_results.R")
source("../R/plotmeanSE.R")

## ----demotypeM,echo=FALSE,fig.width=7,fig.height=7,include=TRUE----------
set.seed(4321)
d<-15
sd<-100
lown<-power.t.test(d=d,sd=sd,power=.10,type="one.sample",alternative="two.sided",strict=TRUE)$n
highn<-power.t.test(d=d,sd=sd,power=.80,type="one.sample",alternative="two.sided",strict=TRUE)$n
nsim<-50
tlow<-thigh<-meanslow<-meanshigh<-CIuplow<-CIlwlow<-CIuphigh<-CIlwhigh<-NULL
critlow<-abs(qt(0.025,df=lown-1))
crithigh<-abs(qt(0.025,df=highn-1))

for(i in 1:nsim){
  x<-rnorm(lown,mean=d,sd=sd)
  meanslow[i]<-mean(x)
  tlow[i]<-t.test(x)$statistic
  CIuplow[i]<-mean(x)+critlow*sd(x)/sqrt(length(x))
  CIlwlow[i]<-mean(x)-critlow*sd(x)/sqrt(length(x))
  x<-rnorm(highn,mean=d,sd=sd)
  meanshigh[i]<-mean(x)
  thigh[i]<-t.test(x)$statistic
  CIuphigh[i]<-mean(x)+crithigh*sd(x)/sqrt(length(x))
  CIlwhigh[i]<-mean(x)-crithigh*sd(x)/sqrt(length(x))
}

 
siglow<-ifelse(abs(tlow)>abs(critlow),"p<0.05","p>0.05")
sighigh<-ifelse(abs(thigh)>abs(crithigh),"p<0.05","p>0.05")

summarylow<-data.frame(means=meanslow,significance=siglow, CIupper=CIuplow, CIlower=CIlwlow)
summaryhigh<-data.frame(index=1:nsim,means=meanshigh,significance=sighigh, CIupper=CIuphigh, CIlower=CIlwhigh)


# re-order data by mean effect size
summarylow<-summarylow[order(summarylow$means), ]
summarylow$index<-1:nrow(summarylow)
summaryhigh<-summaryhigh[order(summaryhigh$means), ]
summaryhigh$index<-1:nrow(summaryhigh)

p_low<-ggplot(summarylow, aes(y=means, x=index,
                              shape=significance,  ymax=CIupper, ymin=CIlower)) + 
  geom_pointrange()+
#  coord_flip()+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(1, 17))+
  magnifytext(sze=22)+ 
  geom_hline(yintercept=15) +
  theme_bw() + 
    scale_x_continuous(name = "")+
  scale_y_continuous(name = "means",limits=c(-100,110))+
  labs(title="Effect 15 ms, SD 100, \n n=20, power=0.10")+
  theme(legend.position="none")+geom_hline(yintercept=0, linetype="dotted")

p_hi<-ggplot(summaryhigh, aes(y=means, x=index,
                              shape=significance, ymax=CIupper, ymin=CIlower)) + 
  geom_pointrange()+
#  coord_flip()+
  geom_point(size=2.5)+
  scale_shape_manual(values=c(1, 17))+
    scale_x_continuous(name = "Sample id")+
  magnifytext(sze=22)+ 
  geom_hline(yintercept=15) +
  theme_bw() + 
  scale_y_continuous(name = "means",limits=c(-100,110))+
  labs(title="Effect 15 ms, SD 100, \n n=350, power=0.80")+
  theme(legend.position=c(0.8,0.3))+geom_hline(yintercept=0, linetype="dotted")

multiplot(p_low,p_hi,cols=1)

## ----LKpredictions,cache=FALSE,echo=FALSE,include=TRUE,fig.width=8,fig.height=4----
LKsurprisal <- data.frame(
  pred = c(700,600,600,500),
  cond = c("a", "b", "c", "d")
)

LKmemory <- data.frame(
  pred = c(500,600,600,700),
  cond = c("a", "b", "c", "d")
)

plot1<-plotpredictions(dat=LKsurprisal,maintitle="Predictions of \n the expectation-based account")
plot2<-plotpredictions(dat=LKmemory,maintitle="Predictions of \n the memory account",ylabel="")
multiplot(plot1,plot2,cols=2)

## ----lkjvisual,cache=TRUE,echo=FALSE,include=FALSE-----------------------
fake_data <- list(x = rnorm(30000,0,1),N = 30000, R = 2) 

stancode <- "
data {
  int<lower=0> N; 
  real x[N]; 
  int R;
  }
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  x ~ normal(mu,sigma);  
}
generated quantities {
  corr_matrix[R] LKJ05;
  corr_matrix[R] LKJ1;
  corr_matrix[R] LKJ2;
  corr_matrix[R] LKJ4;
  LKJ05 = lkj_corr_rng(R,.1);
  LKJ1 = lkj_corr_rng(R,.5);
  LKJ2 = lkj_corr_rng(R,1);
  LKJ4 = lkj_corr_rng(R,2);
}
"

fitfake <- stan(model_code = stancode, pars = c("LKJ05","LKJ1","LKJ2","LKJ4"),
                data = fake_data, chains = 4, 
                iter = 2000)

corrs<-extract(fitfake,pars=c("LKJ05[1,2]","LKJ1[1,2]","LKJ2[1,2]","LKJ4[1,2]"))
corrs<-data.frame(corrs)
colnames(corrs)<-c("lkj05","lkj1","lkj2","lkj4")

## ----figpriors,echo=FALSE,include=TRUE,warning=FALSE,message=FALSE,fig.width=6,fig.height=4----

lkjplot1<-ggplot(corrs, aes(lkj05)) +
  geom_density(adjust=2)+xlab(expression(rho))+ggtitle("nu=0.1")+theme_bw()
lkjplot2<-ggplot(corrs, aes(lkj1)) +
  geom_density(adjust=3)+xlab(expression(rho))+ggtitle(expression(nu==0.5))+theme_bw()+ylim(0, 0.6)
lkjplot3<-ggplot(corrs, aes(lkj2)) +
  geom_density(adjust=3)+xlab(expression(rho))+ggtitle("nu=1")+theme_bw()
lkjplot4<-ggplot(corrs, aes(lkj2)) +
  geom_density(adjust=3)+xlab(expression(rho))+ggtitle(expression(nu==2))+theme_bw()+ylim(0, 0.6)
multiplot(lkjplot2,lkjplot4,cols=2)

## ----nonsigETdepmeasuresLKE2Replication,include=FALSE,eval=FALSE,echo=FALSE,message=FALSE,warning=FALSE----
## # Eyetracking Replication of LK Experiment 2
## # regression probability, skipping probability
## 
## # data ET Replication of LK Expt 2, reload to create a standalone chunk:
## dat<-read.table("../data/E4ETlevykellerExp2.txt",header=T)
## # subset data
## dat<-subset(dat,condition!="f" & condition!="p")
## 
## 
## # add column skipping probability (1 if word skipped, 0 otherwise)
## dat$SKP <- ifelse(dat$TFT>0,0,1)
## # add column regression probability (RBRC is no of regression from word n before continuing to the right.
## # 1 if regression occured, 0 otherwise.)
## dat$FPRP <- as.logical(dat$RBRC)
## 
## # contrast coding
## dat$dat<-ifelse(dat$condition%in%c("a","b"),1/2,-1/2)
## dat$adj<-ifelse(dat$condition%in%c("b","d"),-1/2,1/2)
## dat$int<-ifelse(dat$condition%in%c("b","c"),-1/2,1/2)
## 
## dat$roi<-factor(dat$roi)
## 
## region<-ifelse(dat$roi==23,"npacc",
##                ifelse(dat$roi==25,"verb",
##                       ifelse(dat$roi==27,"verb1","noncritical")))
## 
## dat$region<-region
## dat<-subset(dat,region!="noncritical")
## dat$region<-factor(dat$region)
## 
## # subset critical and postcritical region
## verb <-(subset(dat,region=="verb"))
## verb1 <-(subset(dat,region=="verb1"))
## 
## ## FPRP (first-pass regression prob) crit
## #summary(mFPRP <- glmer(FPRP~dat+adj+int+(1|subject)+(1|itemid), verb, family=binomial))
## 
## ## Skip prob crit
## #summary(mSKP <- glmer(SKP~dat+adj+int+(1|subject)+(1|itemid), verb, family=binomial))
## 
## 
## # FPRP crit
## subj <- as.integer(factor(verb$subject))
## N_subj <- length(unique(subj))
## item <- as.integer(factor(verb$itemid))
## N_items <- length(unique(item))
## 
## 
## X <- unname(model.matrix(~ 1 + dat + adj + int, verb))
## attr(X, which="assign") <- NULL
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb$FPRP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mFPRPStan <- stan(file = "StanModels/logitmaximal.stan",
##                data = stanDat,
##                iter = 2000,
##                chains = 4)
## 
## 
## # FPRP postcrit
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb1$FPRP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mFPRPpost <- stan(file = "StanModels/logitmaximal.stan",
##                   data = stanDat,
##                   iter = 2000,
##                   chains = 4)
## 
## 
## 
## # Skipping prob crit
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb$SKP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mSKPStan <- stan(file = "StanModels/logitmaximal.stan",
##                   data = stanDat,
##                   iter = 2000,
##                   chains = 4)
## 
## 
## 
## 
## # Skipping prob postcrit
## 
## # 2. Make Stan data (list)
## stanDat <- list(N = nrow(X),
##                 P = ncol(X),
##                 n_u = ncol(X),
##                 n_w = ncol(X),
##                 X = X,
##                 Z_u = X,
##                 Z_w = X,
##                 J = N_subj,
##                 K = N_items,
##                 prob = verb1$SKP,
##                 subj = subj,
##                 item = item)
## 
## 
## # 3. Fit the model.
## 
## mSKPpost <- stan(file = "StanModels/logitmaximal.stan",
##                  data = stanDat,
##                  iter = 2000,
##                  chains = 4,
##                  control = list(adapt_delta=0.99))
## 
## 

## ----originalLKdataExp1, include=FALSE, cache=FALSE,echo=FALSE,warning=FALSE----

##### ORIGINAL LEVY KELLER DATA EXPERIMENT 1 #####
# CRITICAL REGION (REGION 7)  #

reading_time <- read.table('../data/exp1_tt_r.res', header=TRUE)
#head(reading_time)

condition<-ifelse(reading_time$dat=="sub" & reading_time$adj=="sub","a",
                  ifelse(reading_time$dat=="sub" & reading_time$adj=="main","b",
                         ifelse(reading_time$dat=="main" & reading_time$adj=="sub","c", 
                                ifelse(reading_time$dat=="main" & reading_time$adj=="main","d","NA"))))

reading_time$condition<-factor(condition)


# contrast coding: 
reading_time$dat<-ifelse(reading_time$condition%in%c("a","b"),1/2,-1/2)
reading_time$adj<-ifelse(reading_time$condition%in%c("b","d"),-1/2,1/2)
reading_time$int<-ifelse(reading_time$condition%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5    
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5

# Thus, ME positive coefficient= longer RTs/slowdown when DAT bzw. PP ADJ in subordinate clause; 
# negative coefficient= faster RTs/speed-up when DAT bzw. PP ADJ in subordinate clause. 
# Interaction positve coefficient= longer RTs/slowdown when DAT/PP ADJ are in the same - subordinate OR main - clause. 

# remove zeros
reading_time_nozeros <- reading_time[reading_time$region7 != 0,]


# model (log tft) at the critical region (orgininally raw rts, see residuals)
#mLKE1crit  <- lmer(log(region7) ~ dat*adj + (dat*adj|subj) + (dat*adj|item), data=reading_time_nozeros)
#summary(mLKE1crit)

# I subset the data frame here (LK1 conds c and d only, such that I can bind it with conditions c and d from LK2).

dataLK1cd<-subset(reading_time, condition!="a" & condition!="b")

## ----AnalysisLK1critical,results='hide',include=FALSE,cache=TRUE,echo=FALSE,warning=FALSE----


# ANALYSIS of TRT original LK1 
# CRITICAL REGION (main verb)
stanDatLKE1<-createStanDat(d=reading_time_nozeros,
              rt=reading_time_nozeros$region7,
              form=as.formula("~1+dat+adj+int"))

LKE1 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE1,
                    iter = 2000, 
                    chains = 4)

LKE1_res<-stan_results(m=LKE1,params=c("Dat","Adj","DatxAdj"))

## original:
#interact  <- lmer(region7 ~ dat+adj+int + (dat+adj+int|subj) + (dat+adj+int|item), data=reading_time_nozeros)
#summary(interact)

## ----AnalysisLK1postcritical,include=FALSE,cache=TRUE,echo=FALSE,warning=FALSE----
# ANALYSIS of TRT original LK1 
# POSTCRITICAL REGION

reading_time_nozeros <- reading_time[reading_time$region8 != 0,]

#mLKE1post  <- lmer(log(region8) ~ dat*adj + (dat+adj|subj) + (dat+adj|item), data=reading_time_nozeros)
#model changed to above (failed to converge before (the original model))

stanDatLKE1post<-createStanDat(d=reading_time_nozeros,
              rt=reading_time_nozeros$region8,
              form=as.formula("~1+dat+adj+int"))

LKE1post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE1post,
                    iter = 2000, 
                    chains = 4)

LKE1post_res<-stan_results(m=LKE1post,params=c("Dat","Adj","DatxAdj"))

## ----originalLKdataExp2,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE----

##### ORIGINAL LEVY KELLER DATA EXPERIMENT 2 #####

# CRITICAL REGION  (REGION 8) #

reading_time <- read.table('../data/exp3_tt_r.res', header=TRUE)


condition<-ifelse(reading_time$dat=="sub" & reading_time$adj=="sub","a",
                  ifelse(reading_time$dat=="sub" & reading_time$adj=="main","b",
                         ifelse(reading_time$dat=="main" & reading_time$adj=="sub","c", 
                                ifelse(reading_time$dat=="main" & reading_time$adj=="main","d","NA"))))

reading_time$condition<-factor(condition)


# contrast coding: 
reading_time$dat<-ifelse(reading_time$condition%in%c("a","b"),1/2,-1/2)
reading_time$adj<-ifelse(reading_time$condition%in%c("b","d"),-1/2,1/2)
reading_time$int<-ifelse(reading_time$condition%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5


# exlude zeros from tfts
reading_time_nozeros <- reading_time[reading_time$region8 != 0,]

#mLKE2crit  <- lmer(log(region8) ~ dat*adj + (dat*adj|subj) + (dat*adj|item), data=reading_time_nozeros)
#summary(mLKE2crit)

# Subset the data frame (conds c and d only such that one can bind it with conditions c and d from LK1).
# LK2 c, d do not need to be relabelled

dataLK2cd<-subset(reading_time, condition!="a" & condition!="b")

## ----AnalysisLK2critical,results='hide',include=FALSE, cache=TRUE,echo=FALSE,warning=FALSE----

# ANALYSIS of TRT original LK2
# CRITICAL REGION (verb+aux)

stanDatLKE2<-createStanDat(d=reading_time_nozeros,
                       rt=reading_time_nozeros$region8,
                       form=as.formula("~1+dat+adj+int"))

#str(stanDat)
LKE2 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE2,
                    iter = 2000, 
                    chains = 4)

LKE2_res<-stan_results(m=LKE2,params=c("Dat","Adj","DatxAdj"))

## ----AnalysisLK2postcritical,include=FALSE,results='hide',cache=TRUE,echo=FALSE,warning=FALSE----

# EXP 2 REGION 9 POSTCRITICAL REGION #

reading_time_nozeros <- reading_time[reading_time$region9 != 0,]

#mLKE2post  <- lmer(log(region9) ~ dat*adj + (dat*adj|subj) + (dat+adj|item), data=reading_time_nozeros)
#summary(mLKE2post)


stanDatLKE2post<-createStanDat(d=reading_time_nozeros,
                       rt=reading_time_nozeros$region9,
                       form=as.formula("~1+dat+adj+int"))

LKE2post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatLKE2post,
                    iter = 2000, 
                    chains = 4)

LKE2post_res<-round(stan_results(m=LKE2post,
                           params=c("Dat","Adj","DatxAdj")))

## ----ResponseAccSPRLK1data,include=FALSE,echo=FALSE,warning=FALSE,message=FALSE----
dat<-read.table("../data/E1SPRlevykellerExp1.txt",header=T)
datModel = read.csv("../../recursive-prd/output/E1SPRlevykellerExp1.txt_7092426", sep="\t")

dat$LineNumber = (1:nrow(dat)) - 1
dat = merge(dat, datModel, by=c("LineNumber"))
dat = dat[!grepl("OOV", dat$RegionLSTM),]
dat = dat[dat$RegionLSTM!="",]
dat = dat[dat$subj == 1,]
#head(dat)
#str(dat)

dat$item<-factor(dat$item)
dat$subj<-factor(dat$subj)

# subset experimental items (exlude filler, practice items)
E1spr<-subset(dat,roi!="?" & expt=="LKrep")

# in original LK1, postcrit region is "und so/und damit", hence, here sum rts of verb1+verb2

E1spr$region<-ifelse(E1spr$roi==13,"verb",
                     ifelse(E1spr$roi==14, "verb1", 
                            ifelse(E1spr$roi==15, "verb2", "noncritical")))


E1sprCRIT<-subset(E1spr, region=="verb")
E1sprPOST1<-subset(E1spr, region=="verb1")
E1sprPOST2<-subset(E1spr, region=="verb2")

#E1sprCRIT<-E1sprCRIT[,c(1,3,4,7,8)]
#E1sprPOST1<-E1sprPOST1[,c(1,3,4,7,8)]
#E1sprPOST2<-E1sprPOST2[,c(1,3,4,7,8)]


E1sprPOST1$rt<-E1sprPOST1$rt + E1sprPOST2$rt
E1sprPOST1$Surprisal<-E1sprPOST1$Surprisal + E1sprPOST2$Surprisal

# inspect data frame/double-check the right rts were added
#E1sprPOST<-cbind(E1sprPOST1,E1sprPOST2)
#E1sprPOST$rtPOST<-E1sprPOST1$rt + E1sprPOST2$rt

E1spr<-rbind(E1sprCRIT,E1sprPOST1)

# contrast coding: 
E1spr$dat<-ifelse(E1spr$cond%in%c("a","b"),1/2,-1/2)
E1spr$adj<-ifelse(E1spr$cond%in%c("b","d"),-1/2,1/2)
E1spr$int<-ifelse(E1spr$cond%in%c("b","c"),-1/2,1/2)

                 ## ME DAT ## ME PP-ADJ ## INT
# a DAT-SC; PP-SC    0.5         0.5       0.5
# b DAT-SC; PP-MC    0.5        -0.5      -0.5
# c DAT-MC; PP-SC   -0.5         0.5      -0.5
# d DAT-MC; PP-MC   -0.5        -0.5       0.5


# subset critical, postcritical 
verb<-subset(E1spr,region=="verb")
verb1<-subset(E1spr,region=="verb1")

## ----AnalysisSPRLK1critical,cache=TRUE,include=FALSE---------------------

# E1 SPR rt at critical 
#head(verb)
stanDatSPRLK1<-createStanDat(d=verb,rt=verb$rt,
                             form=as.formula("~1+dat+adj+int"))


summary(lm(Surprisal ~ 1+dat+adj+int, data=verb))
summary(lm(Surprisal ~ 1+dat+adj+int, data=verb1))

SPRLK1 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRLK1,
                    iter = 2000, 
                    chains = 4)

# check
#m1 <- lmer(log(rt)~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb)
#summary(m1)

stanDatSPRE1post<-createStanDat(d=verb1,rt=verb1$rt,form=as.formula("~1+dat+adj+int"))

SPRE1post <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRE1post,
                    iter = 2000, 
                    chains = 4)

# check
#m1post <- lmer(log(rt)~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb1)
#summary(m1post)

SPRE1post_res<-round(stan_results(m=SPRE1post,params=pars))

## ----AnalysisSPRLK1postcriticalraw,echo=FALSE----------------------------
#SPRLK1postcritm1post <- lmer(rt~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb1)

## ----AccE3data,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE--------
dat<-read.table("../data/E3SPRlevykellerExp2.txt",header=T)
datModel = read.csv("../../recursive-prd/output/E3SPRlevykellerExp2.txt_7092426", sep="\t")

dat$LineNumber = (1:nrow(dat)) - 1
dat = merge(dat, datModel, by=c("LineNumber"))
dat = dat[!grepl("OOV", dat$RegionLSTM),]
dat = dat[dat$RegionLSTM!="",]
dat = dat[dat$subj == 1,]

dat$item<-factor(dat$item)
dat$subj<-factor(dat$subj)

E3spr<-subset(dat,roi!="?" & expt=="LKrep")


# contrast coding (same as E1 above)
E3spr$dat<-ifelse(E3spr$cond%in%c("a","b"),1/2,-1/2)
E3spr$adj<-ifelse(E3spr$cond%in%c("b","d"),-1/2,1/2)
E3spr$int<-ifelse(E3spr$cond%in%c("b","c"),-1/2,1/2)

## subset critical region, postcritical region 
verb<-subset(E3spr,region=="verb")
verb1<-subset(E3spr,region=="verb1")
#head(verb)


summary(lm(Surprisal~dat+adj+int, data=verb))
summary(lm(Surprisal~dat+adj+int, data=verb1))

# mean RTs at critical region:
#with(verb,round(tapply(rt,cond,mean)))
# mean RTs at postcritical region:
#with(verb1,round(tapply(rt,cond,mean)))


## ----AnalysisRTcritE3,cache=TRUE,echo=FALSE,include=FALSE----------------

stanDatSPRE2<-createStanDat(d=verb,rt=verb$rt,
              form=as.formula("~1+dat+adj+int"))

SPRE2 <- stan(file = "StanModels/maxModel.stan", 
                    data = stanDatSPRE2,
                    iter = 2000, 
                    chains = 4)

# check
#m2 <- lmer(log(rt)~dat+adj+int+(1+dat+adj+int||subj)+(1+dat+adj+int||item),verb)
#summary(m2)

SPRE2_res<-round(stan_results(m=SPRE2,params=pars))

## ----E5ResponseAccdata,include=FALSE,cache=FALSE,echo=FALSE,warning=FALSE----
dat<-read.table("../data/E5SPRlevykellerExp12.txt",header=T)
datModel = read.csv("../../recursive-prd/output/E5SPRlevykellerExp12.txt_7092426", sep="\t")

dat$LineNumber = (1:nrow(dat)) - 1
dat = merge(dat, datModel, by=c("LineNumber"))
dat = dat[!grepl("OOV", dat$RegionLSTM),]
dat = dat[dat$RegionLSTM!="",]
dat = dat[dat$subj == 1,]

dat$item<-factor(dat$item)
dat$subj<-factor(dat$subj)


E5spr<-subset(dat,roi!="?" & expt=="LKrep")

#round(with(subset(E5spr,region=="verb1"),tapply(rt,cond,mean)))
#round(with(subset(E5spr,region=="verb2"),tapply(rt,cond,mean)))


# ROIs: crit: "verb" (a,b = "versteckt, c,d = "versteckt hat"), 
# postcrit: "verb1" (a,b = "und damit", c,d = "den Besuch")
# postcrit for a,b not presented together in Linger, therefore, merge here (verb1+verb2)

#datab<-subset(E5spr,cond!="c" & cond!="d")

#E5sprCRIT<-subset(datab,region=="verb")
#E5sprPOST1<-subset(datab,region=="verb1")

#E5sprPOST2<-subset(datab,region=="verb2")

#E5sprCRIT<-E5sprCRIT[,c(1,3,4,7,8)]
#E5sprPOST1<-E5sprPOST1[,c(1,3,4,7,8)]
#E5sprPOST2<-E5sprPOST2[,c(1,3,4,7,8)]

#E5sprPOST1$rt<-E5sprPOST1$rt + E5sprPOST2$rt

#datab<-rbind(E5sprCRIT,E5sprPOST1)

#datcd<-subset(E5spr,cond!="a" & cond!="b")
#datcd<-subset(datcd,region==c("verb","verb1"))   # "verb1 region",e.g. "den_Konkurs" was shown as one roi in Linger
#datcd<-datcd[,c(1,3,4,7,8)]

#E5spr<-rbind(datab,datcd)

#boxplot(log(rt)~cond,subset(E5spr,region=="verb1"))


# contrast coding
E5spr$load<-ifelse(E5spr$cond%in%c("a","b"),-1/2,1/2)
E5spr$dist<-ifelse(E5spr$cond%in%c("a","c"),-1/2,1/2)
E5spr$int<-ifelse(E5spr$cond%in%c("a","d"),-1/2,1/2)

                         ## ME LOAD ## ME DIST  ## INT
# a DAT-MC;     PP-SC      -0.5        -0.5      -0.5   (originally E1 LK13 cond c)
# b DAT-MC;     PP-MC      -0.5         0.5       0.5   (orininally E1 LK13 cond d)
# c DAT-MC emb; PP-SC emb   0.5        -0.5       0.5   (originally E2 LK13 cond c)
# d DAT-MC emb; PP-MC emb   0.5         0.5      -0.5   (originally E2 LK13 cond d)



# subset crit, postcrit
verb<-subset(E5spr,region=="verb")
verb1<-subset(E5spr,region=="verb1")




summary(lm(Surprisal ~ 1+load+dist+int, data=verb))
summary(lm(Surprisal ~ 1+load+dist+int, data=verb1))

## ----sprmergedraw,echo=FALSE---------------------------------------------
#msprmergedcritraw<-lmer(rt~load + dist + int+(1+load + dist + int||subj)+(1+load + dist + int||item),verb)
#msprmergedpostcritraw<-lmer(rt~load + dist + int+(1+load + dist + int||subj)+(1+load + dist + int||item),verb1)


