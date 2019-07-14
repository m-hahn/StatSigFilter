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

## ----ResponseAccSPRLK1data,include=FALSE,echo=FALSE,warning=FALSE,message=FALSE----
dat<-read.table("../data/E1SPRlevykellerExp1.txt",header=T)
#datModel = read.csv("../../recursive-prd/output/E1SPRlevykellerExp1.txt_7092426", sep="\t")
#datModel = read.csv("../../recursive-prd/output/E1SPRlevykellerExp1.txt_wiki-german-nospaces-bptt-WHITESPACE-39149757", sep="\t")
datModel = read.csv("../../recursive-prd/output/E1SPRlevykellerExp1.txt_657350374", sep="\t")

dat$LineNumber = (1:nrow(dat)) - 1
dat = merge(dat, datModel, by=c("LineNumber"))
dat = dat[!grepl("OOV", dat$RegionLSTM),]
dat = dat[dat$RegionLSTM!="",]
#dat = dat[dat$subj == 1,]
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


summary(lmer(Surprisal ~ 1+dat+adj+int+(1+dat+adj+int|item), data=verb))
summary(lmer(Surprisal ~ 1+dat+adj+int+(1+dat+adj+int|item), data=verb1)) 

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
#dat = dat[dat$subj == 1,]

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


summary(lmer(Surprisal~dat+adj+int+(1+dat+adj+int|item), data=verb)) # evidence for effects of dat and int
summary(lmer(Surprisal~dat+adj+int+(1+dat+adj+int|item), data=verb1))

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
#datModel = read.csv("../../recursive-prd/output/E5SPRlevykellerExp12.txt_7092426", sep="\t")

models = c(
710899603,
553555034,
657350374,
887261094,
459934402)

datModel = data.frame()
for(model in models) {
    datModel2 = read.csv(paste("../../recursive-prd/output/E5SPRlevykellerExp12.txt_", model, sep=""), sep="\t")
    datModel2$Model = model
    datModel = rbind(datModel, datModel2)
}
dat$LineNumber = (1:nrow(dat)) - 1
dat = merge(dat, datModel, by=c("LineNumber"))
dat = dat[!grepl("OOV", dat$RegionLSTM),]
dat = dat[dat$RegionLSTM!="",]
#dat = dat[dat$subj == 1,]

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

summary(lmer(Surprisal ~ 1+load+dist+int+(1+load+dist+int|item)+(1+load+dist+int|Model), data=verb)) # small evidence for effect of distance (also in the character-based model)

library(brms)

model = brm(Surprisal ~ 1+load+dist+int+(1+load+dist+int|item)+(1+load+dist+int|Model), data=verb)

summary(model)



summary(lmer(Surprisal ~ 1+load+dist+int+(1+load+dist+int|item), data=verb)) # small evidence for effect of distance (also in the character-based model)
summary(lmer(Surprisal ~ 1+load+dist+int+(1+load+dist+int|item), data=verb1))

## ----sprmergedraw,echo=FALSE---------------------------------------------
#msprmergedcritraw<-lmer(rt~load + dist + int+(1+load + dist + int||subj)+(1+load + dist + int||item),verb)
#msprmergedpostcritraw<-lmer(rt~load + dist + int+(1+load + dist + int||subj)+(1+load + dist + int||item),verb1)


