
setwd("~/git/tfg_assess2020")
library(TMB)
library(tidyverse)
library(glmmTMB)
library(frasyr)
# install.packages("DHARMa")
library(DHARMa)
library(mgcv)

source("surfnet_peak_modeling.R")

TMB::compile("surfnet_peak_model.cpp")
dyn.load(dynlib("surfnet_peak_model"))

# data の読み込み ----
# dat0 = read.csv("C:/Users/00007802/Dropbox/tfg/index_age0/surfnet_survey2020.csv", header=T)
# dat0 = read.csv("data/surfnet_survey2020.csv", header=T)
Dat0 = read.csv("data/surfnet_survey2021.csv", header=T)

Dat = Dat0 %>% 
  # mutate(Time = scale(Time)) %>%
  select(Catch,Effort,Year,Time) %>%
  arrange(Year,Time)


## 網の大きさに応じた重みづけの試算 ----

Dat= Dat0 %>% 
  mutate(area_w = case_when(Year %in% c(2004,2005) ~ 5*1.5/4,
                            Year %in% c(2006,2007) ~ 10/4,
                            TRUE ~ 1),
         length_w = case_when(Year %in% c(2004,2005) ~ 5/4,
                              Year %in% c(2006,2007) ~ 10/4,
                              TRUE ~ 1)) %>%
  mutate(Effort = Effort*area_w)   # 面積で重みづけする場合
# mutate(Effort = Effort*length_w)   # 長さで重みづけする場合

Dat1 = Dat %>% 
  # mutate(Time = scale(Time)) %>%
  select(Catch,Effort,Year,Time) %>%
  arrange(Year,Time)
nrow(Dat) #データ数=105

Dat2= Dat0 %>% 
  mutate(area_w = case_when(Year %in% c(2004,2005) ~ 5*1.5/4,
                            Year %in% c(2006,2007) ~ 10/4,
                            TRUE ~ 1),
         length_w = case_when(Year %in% c(2004,2005) ~ 5/4,
                              Year %in% c(2006,2007) ~ 10/4,
                              TRUE ~ 1)) %>%
  # mutate(Effort = Effort*area_w)   # 面積で重みづけする場合
  mutate(Effort = Effort*length_w)   # 長さで重みづけする場合

Dat2 = Dat2 %>% 
  # mutate(Time = scale(Time)) %>%
  select(Catch,Effort,Year,Time) %>%
  arrange(Year,Time)

tapply(Dat$Catch,Dat$Year,length) %>% summary
tapply(Dat$Catch,Dat$Year,length) %>% sd

Dat$Effort %>% summary()
Dat$Effort %>% sd()

knitr::kable(head(Dat))
knitr::kable(head(Dat1))
knitr::kable(head(Dat2))

# Effortが異なる
# rm(dat)
basemodel = spm(dat=Dat,
                family="nbinom2",
                a_key=2,b_key=1,c_key=4,
                silent = FALSE)

areamodel = spm(dat=Dat1,
                family="nbinom2",
                a_key=2,b_key=1,c_key=4,
                silent = FALSE)

lengthmodel = spm(dat=Dat2,
                  family="nbinom2",
                  a_key=2,b_key=1,c_key=4,
                  silent = FALSE)

c(basemodel$BIC,areamodel$BIC,lengthmodel$BIC)
# length modelが少しだけ良いのでDat2を使う


## 網羅的なモデル選択 ----

basemodel = spm(dat=Dat2,
                family="nbinom2",
                a_key=1,b_key=1,c_key=2,
                silent = TRUE)

load("bestmodel.rda")
# bestmodel$input$dat
(g_fit = plot_spm_fit(basemodel,CI=0.9))

basemodel2 = spm(dat=Dat2,
                family="nbinom2",
                a_key=2,b_key=1,c_key=4,
                silent = TRUE)

(g_fit2 = plot_spm_fit(basemodel2,CI=0.8))
#2006年が広い

(g_trend = plot_spm_trend(basemodel))
(g_trend2 = plot_spm_trend(basemodel2))
  

# nbinom2だけで
model_list = expand.grid(family = 2, # 0: poisson, 1: nbinom1, 2: nbinom2
                         a_key = 0:4,
                         b_key = 0:4,
                         c_key = 0:4)

nrow(model_list) #N = 125

res_list=list()
summary_table = c()

minAICc = 9999

# i=32

# i = 1

# i = 28

for (i in 1:nrow(model_list)) {
  input_tmp = basemodel$input
  # basemodel$input$dat$Effort
  # input_tmp$p0_list <- basemodel2$par_list

  if (model_list[i,1] == 0) {
    input_tmp$family = "poisson"
  } else {
    if (model_list[i,1] == 1) {
      input_tmp$family = "nbinom1"
    } else {
      input_tmp$family = "nbinom2"
    }
  }
  input_tmp$a_key = model_list[i,2]
  input_tmp$b_key = model_list[i,3]
  input_tmp$c_key = model_list[i,4]
  res = try(do.call(spm,input_tmp))
  # res$rep
  # res$opt$convergence
  res_list[[i]] <- res
  
  if (class(res) == "try-error") {
    res_vec = rep(NA,8) 
  } else {
    res_vec = c(convergence = res$convergence,
                pdHess = res$pdHess,
                N = res$N,
                npar = res$npar,
                loglik = res$loglik,
                AIC = res$AIC,
                AICc = res$AICc,
                BIC = res$BIC) 
    if (res$convergence==0 && res$pdHess) {
      if (res$AICc < minAICc) minAICc = res$AICc
    }
  }
  summary_table = cbind(summary_table,c(id = i, unlist(model_list[i,]),res_vec))
  
  cat("[",i,"] ", "pdHess: ",res_vec[2],", AICc: ", res_vec[7]," (min=",minAICc,")\n",sep="")
}


SD_est = sapply(1:ncol(summary_table), function(i) {
  if(class(res_list[[i]])!="try-error") {
    sum(is.nan(res_list[[i]]$rep$sd))==0
  }else{
    FALSE
  } 
})

summary_table2 = t(summary_table) %>% as_tibble() %>%
  mutate(SD_est = SD_est) %>% 
  arrange(convergence,desc(pdHess),desc(SD_est),AIC)

write.csv(summary_table2, file = "summary_AICtable.csv")

summary_table2 = t(summary_table) %>% as_tibble() %>% 
  mutate(SD_est = SD_est) %>% 
  arrange(convergence,desc(pdHess),desc(SD_est),AIC) %>%
  mutate(dAIC = AIC-AIC[1]) %>%
  arrange(convergence,desc(pdHess),desc(SD_est),AICc) %>%
  mutate(dAICc = AICc-AICc[1]) %>% 
  arrange(convergence,desc(pdHess),desc(SD_est),BIC) %>%
  mutate(dBIC = BIC-BIC[1])

write.csv(summary_table2, file = "summary_BICtable.csv")

summary_table3 = summary_table2 %>%
  arrange(convergence,desc(pdHess),desc(SD_est),AIC)

write.csv(summary_table3, file = "summary_AICtable.csv")

summary_table4 = summary_table2 %>%
  arrange(convergence,desc(pdHess),desc(SD_est),AICc)

write.csv(summary_table4, file = "summary_AICctable.csv")

head(summary_table2) #BIC (2,1,4)
head(summary_table3) #AIC (1,2,0)
head(summary_table4) #AICc (2,1,4)

# summary_table3

# AICcとBICは一致する

best_BIC = res_list[[summary_table2$id[1]]]
c(best_BIC$input$a_key,best_BIC$input$b_key,best_BIC$input$c_key)
best_BIC$loglik
best_BIC$rep

best_AIC = res_list[[summary_table3$id[1]]]
c(best_AIC$input$a_key,best_AIC$input$b_key,best_AIC$input$c_key)
best_AIC$loglik
best_AIC$rep

# res_list[[summary_table3$id[2]]]$loglik

best_AIC$convergence
best_AIC$pdHess

sum(is.nan(best_AIC$rep$sd))
sum(is.nan(best_BIC$rep$sd))

g_aic_fit = plot_spm_fit(best_AIC)
g_aic_fit

# model_aic = spm(dat=Dat2,
#                  family="nbinom2",
#                  a_key=0,b_key=0,c_key=1,
#                  silent = TRUE
#                 # ,p0_list=basemodel2$par_list
#                 )
# model_aic$rep
# model_aic$convergence
# model_aic$pdHess
# 
# model_aic$loglik
# best_AIC$loglik
# best_AIC$rep
# 
# model_aic2 = spm(dat=Dat2,
#                 family="nbinom2",
#                 a_key=1,b_key=2,c_key=0,
#                 silent = TRUE,p0_list=basemodel2$par_list)
# model_aic2$rep
# model_aic2$convergence
# model_aic2$pdHess

# model_aic3 = spm(dat=Dat2,
#                  family="nbinom2",
#                  a_key=0,b_key=2,c_key=3,
#                  silent = TRUE
#                  ,p0_list=best_BIC$par_list
#                  )
# 
# model_aic3$rep
# 
# c(model_aic$AIC,model_aic2$AIC)


best_AICc = res_list[[summary_table4$id[1]]]
c(best_AICc$input$a_key,best_AICc$input$b_key,best_AICc$input$c_key)
best_AICc$loglik
best_AICc$rep

### output AIC & BIC table ----

summary_table = read.csv("summary_AICtable.csv",header=T)
summary_table %>% colnames()
summary_table$convergence %>% unique()
summary_table$SD_est %>% unique()
summary_table2 = summary_table %>% filter(convergence==0 & pdHess==1 & SD_est==TRUE) %>%
  mutate(a = case_when(a_key==0 ~ "FE",
                       a_key==1 ~ "Const",
                       a_key==2 ~ "WN",
                       a_key==3 ~ "AR(1)",
                       a_key==4 ~ "RW"),
         b = case_when(b_key==0 ~ "FE",
                       b_key==1 ~ "Const",
                       b_key==2 ~ "WN",
                       b_key==3 ~ "AR(1)",
                       b_key==4 ~ "RW"),
         c = case_when(c_key==0 ~ "FE",
                       c_key==1 ~ "Const",
                       c_key==2 ~ "WN",
                       c_key==3 ~ "AR(1)",
                       c_key==4 ~ "RW"))

aic_table = summary_table2 %>% arrange(AIC) %>%
  mutate(Rank = 1:n()) %>%
  filter(Rank < 6) %>%
  select(Rank,npar,loglik,dAIC,a,b,c) %>% 
  mutate(loglik = round(loglik,2),dAIC=round(dAIC,2)) %>%
  rename(df=npar,"Log-likelihood"=loglik,delta=dAIC)


bic_table = summary_table2 %>% arrange(BIC) %>%
  mutate(Rank = 1:n()) %>%
  filter(Rank < 6) %>%
  select(Rank,npar,loglik,dBIC,a,b,c) %>% 
  mutate(loglik = round(loglik,2),dBIC=round(dBIC,2)) %>%
  rename(df=npar,"Log-likelihood"=loglik,delta=dBIC)

aicc_table = summary_table2 %>% arrange(AICc) %>%
  mutate(Rank = 1:n()) %>%
  filter(Rank < 6) %>%
  select(Rank,npar,loglik,dAICc,a,b,c) %>% 
  mutate(loglik = round(loglik,2),dBIC=round(dAICc,2)) %>%
  rename(df=npar,"Log-likelihood"=loglik,delta=dAICc)


both_table = full_join(aic_table , bic_table)
write.csv(both_table, file="AIC-BIC_table.csv",row.names=FALSE)


# bestmodel=best_AICc

### cross validation ----

# if (0) {
  model_list = list(best_BIC,best_AIC)
  
  # sapply(1:30, function(i) model_list[[i]]$pdHess)
  
  cv_seed1 = spm_CV(model_list,k=nrow(Dat2),seed=1,type=3)
  # cv_seed1 = spm_CV(model_list,k=10,seed=1,type=3)
  # # cv_seed1 = spm_CV(model_list,k=10,seed=1,type=2) #これだとRMSEとMAEでは逆転する
  
  # cv_seed1$summary %>% View()
  cv_summary1 = cv_seed1$summary %>% arrange(mean_rmse)
  cv_summary1 = cv_seed1$summary %>% arrange(mean_mae)
  cv_summary1 = cv_seed1$summary %>% arrange(mean_ls)
  # cv_seed1$all %>% View
  # View(cv_summary1)
  
  # colnames(cv_seed1$all)
  
  cv_summary2 = cv_seed1$all %>% as_tibble() %>% group_by(model) %>%
    summarise(RMSE = sqrt(mean(rmse^2)),MAE = mean(mae), MAD= median(mad), LogScore=mean(ls),SD_LogScore=sd(ls)) %>%
    mutate(Model = if_else(model==1,"Best BIC","Best AIC")) %>%
    select(-model) %>%
    select(Model,everything())
    
  bestmodel = model_list[[cv_summary1$model[1]]]
  
  # View()
  
  # write.csv(cv_summary1, file = "loocv_summary.csv")
  write.csv(cv_summary2, file = "loocv_summary.csv")
  
  write.csv(cv_seed1$all, file = "loocv_all.csv")
  
# }

# save(bestmodel,file="bestmodel.rda")
load("bestmodel.rda")
bestmodel$input$dat %>% nrow

(g_fit_best = plot_spm_fit(bestmodel,CI=0.8,ncol=6))

ggsave(g_fit_best, file = "CPUE_fit_bestmodel_rev.png",dpi=600,unit="mm",height=120,width=240)

Index = bestmodel$index

write.csv(Index,file = "index_bestmodel.csv")

g_trend_best = plot_spm_trend2(bestmodel)
g_trend_best$graph

ggsave(g_trend_best$graph, file = "index_trend.png",dpi=600,unit="mm",height=150,width=240)


### 

### model diagnostics ----

set.seed(1010)
catch_obs = Dat2$Catch
catch_add = catch_obs + runif(length(catch_obs),-0.5,0.5)
sim_catch <- sim_catch_add <- NULL
nsim = 10000
for (i in 1:nsim) {
  sim = bestmodel$obj$simulate()$Catch
  sim_catch = rbind(sim_catch,sim)
  sim_catch_add = rbind(sim_catch_add, sim + runif(length(sim),-0.5,0.5))
  # catch_obs_add = rbind(catch_obs_add, catch_obs + runif(length(catch_obs),-0.5,0.5))
}

# dim(sim_catch)
# ?createDHARMa
tmp = createDHARMa(t(sim_catch),Dat2$Catch,fittedPredictedResponse = bestmodel$Catch_pred,integerResponse = T,seed=1)
plot(tmp)
plot(tmp,quantreg=FALSE)

Resid_test = testResiduals(tmp)
Resid_test$dispersion
Resid_test$outliers
Resid_test$uniformity

ZI_test = testZeroInflation(tmp)
ZI_test$p.value

# ?testDispersion
# Disp_test = test

plot(tmp$scaledResiduals~std_resid )
# dim(catch_obs_add)

# j=1
# sim_catch_add[,j]
# std_resid = sapply(1:length(catch_obs), function(j) sum(sim_catch_add[,j]<catch_obs_add[,j])/nsim)

# hist(std_resid)

# # j=1
# std_resid = sapply(1:length(catch_obs), function(j) sum(sim_catch_add[,j]<catch_add[j])/nsim )
std_resid = tmp$scaledResiduals
zero_count = sapply(1:nsim, function(j) sum(sim_catch[j,]==0))

hist(std_resid)

## 

# nrow(bestmodel$dat_pred)

df = bestmodel$dat_pred %>% mutate(resid = std_resid)

check_ks <- ks.test(std_resid,y="punif",min=0,max=1)
check_ks # Not significant

pvalue = check_ks$p.value

# df = df %>% mutate(label=sprintf("KS test: P == %.2f",pvalue))
df_zero = data.frame(value = zero_count)

hist(zero_count)

# p_ks = data.frame(label = paste0("KS test: P = ", round(pvalue,3)),x=0,y=10)
p_ks = data.frame(label = paste0("KS test: P = ", round(Resid_test$uniformity$p.value,3)),x=0,y=10)

# p_ks = data.frame(label = paste0("KS test: P = ", round(pvalue,3)),x=0,y=15)

g_resid_hist = ggplot(data=df,aes(x=resid)) +
  geom_histogram(breaks=seq(0,1,by=0.05))+
  geom_text(data=p_ks,parse=FALSE,aes(x=x,y=y,label=label,hjust=0,vjust=1),size=4) +
  xlab("Standardized residual")+ylab("Count")+
  ggtitle("(a) Histogram of standardized residuals")+
  frasyr::theme_SH()
g_resid_hist

g_qq = ggplot(data=df) + 
  geom_abline(intercept=0,slope=1,colour="red",size=1)+
  geom_qq(aes(sample=resid),distribution = stats::qunif,
          dparams = list(min=0,max=1),size=2)+
  xlim(0,1)+ylim(0,1)+
  ylab("Standardized residual")+xlab("Theoretical")+
  ggtitle("(b) QQ plot of standardized residuals")+
  theme_SH()
g_qq

colnames(df)
range(std_resid)

df = df %>% mutate(Rank = rank(Catch_pred)) %>%
  mutate(Standardized_Rank = (Rank-min(Rank))/(max(Rank)-min(Rank)))

Rank = function(x) (rank(x)-min(rank(x)))/(max(rank(x))-min(rank(x)))

df = df %>% 
  mutate(Time_rank = Rank(Time),Year_rank = Rank(Year))

# hist((rank(df$Time)-min(rank(df$Time)))/(max(rank(df$Time))-min(rank(df$Time))))

logit=function(x)log(x/(1-x))

m1 = gam(logit(resid)~s(Rank),data=df)
m1 = gam(resid~s(Rank),data=df,family=betar(link="logit"))
# names(summary(m1))
p1 = summary(m1)[["s.table"]][1,4]
p1 = data.frame(label = paste0("P = ", round(p1,3)),x=1,y=1.0)

# tmp2 = loess(logit(resid)~Rank,data=df)
# summary(tmp2)

# summary(glm(logit(resid)~Rank_time,data=df))
# summary(glm(logit(resid)~Time,data=df))
# 
# summary(gam(logit(resid)~s(Rank_time),data=df))
# m2 = gam(logit(resid)~s(Time),data=df)
m2=gam(resid~s(Time),data=df,family=betar(link="logit"))
summary(m2)
# m2 = gam(resid~s(Time),data=df)
# hist(m2$residuals)
# m2 = gam(logit(resid)~s(Time_rank),data=df)
# summary(m2)
p2 = summary(m2)[["s.table"]][1,4]
p2 = data.frame(label = paste0("P = ", round(p2,3)),x=max(df$Time),y=1.0)

m3 = gam(logit(resid)~s(Year),data=df)
# m3 = gam(logit(resid)~s(Year_rank),data=df)
m3 = gam(resid~s(Year),data=df,family=betar(link="logit"))

# plot(predict(m3,type="response")~df$Year,ylim=c(0,1))

p3 = summary(m3)[["s.table"]][1,4]
p3 = data.frame(label = paste0("P = ", round(p3,3)),x=max(df$Year),y=1.0)

g_scatter_smooth = ggplot(data=df,aes(x=Standardized_Rank,y=resid)) +
# g_scatter_smooth = ggplot(data=df,aes(x=Catch_pred,y=log(resid/(1-resid)))) +
  geom_hline(yintercept=0.5,size=1, linetype="dashed",colour="red")+
  geom_point(size=2)+
  # geom_smooth(method="gam",size=1.5,method.args = list(family = gaussian(link = "logit")))+
  geom_smooth(method="gam",size=1.5,method.args = list(family = betar(link="logit")))+
  # geom_smooth(method="loess",size=1.5)+
  xlab("Predicted catch (rank transformed)")+ylab("Standardized residual")+
  ggtitle("(d) Standardized residuals vs. Prediction")+
  theme_SH()+
  ylim(0,1.05)+
  geom_text(data=p1,parse=FALSE,aes(x=x,y=y,label=label,hjust=1,vjust=0),size=4)
g_scatter_smooth

# df = df %>% mutate(Rank_time = rank(Time)) %>%
#   mutate(Standardized_Rank_time = (Rank_time-min(Rank_time))/(max(Rank_time)-min(Rank_time)))

g_scatter_smooth2 = ggplot(data=df,aes(x=Time,y=resid)) +
  # g_scatter_smooth = ggplot(data=df,aes(x=Catch_pred,y=log(resid/(1-resid)))) +
  geom_hline(yintercept=0.5,size=1, linetype="dashed",colour="red")+
  geom_point(size=2)+
  # geom_smooth(method="gam",size=1.5,method.args = list(family = gaussian(link = "logit")))+
  geom_smooth(method="gam",size=1.5,method.args = list(family = betar(link = "logit")))+
  # geom_smooth(method="loess",size=1.5)+
  xlab("Days elapsed since May 1")+ylab("Standardized residual")+
  ggtitle("(f) Standardized residuals vs. Date")+
  theme_SH()+ylim(0,1.05)+
  geom_text(data=p2,parse=FALSE,aes(x=x,y=y,label=label,hjust=1,vjust=0),size=4)
g_scatter_smooth2

g_scatter_smooth3 = ggplot(data=df,aes(x=Year,y=resid)) +
  # g_scatter_smooth = ggplot(data=df,aes(x=Catch_pred,y=log(resid/(1-resid)))) +
  geom_hline(yintercept=0.5,size=1, linetype="dashed",colour="red")+
  geom_point(size=2)+
  # geom_smooth(method="gam",size=1.5,method.args = list(family = gaussian(link = "logit")))+
  geom_smooth(method="gam",size=1.5,method.args = list(family = betar(link = "logit")))+
  # geom_smooth(method="loess",size=1.5)+
  xlab("Year")+ylab("Standardized residual")+
  ggtitle("(e) Standardized residuals vs. Year")+
  theme_SH()+
  ylim(0,1.05)+
  geom_text(data=p3,parse=FALSE,aes(x=x,y=y,label=label,hjust=1,vjust=0),size=4)
g_scatter_smooth3


# pvalue2 = 2*min(sum(df_zero$value<sum(Dat2$Catch==0))/nrow(df_zero),sum(df_zero$value>sum(Dat2$Catch==0))/nrow(df_zero))
pvalue2 = ZI_test$p.value
# hist(df_zero$value)
zero_p = data.frame(label=paste0("P = ",round(pvalue2,3)),x=min(df_zero$value),y=max(table(df_zero)))

g_zero_hist = ggplot(data=df_zero,aes(x=value)) +
  geom_histogram(breaks=min(df_zero$value):max(df_zero$value))+
  geom_vline(xintercept = sum(Dat2$Catch==0),colour="blue",linetype="dashed",size=1.5)+
  geom_text(data=zero_p,parse=FALSE,aes(x=x,y=y,label=label,hjust=0,vjust=1),size=4)+
  xlab("Number of zeros")+ylab("Count")+
  ggtitle("(c) Histogram of zero numbers")+
  theme_SH()
g_zero_hist

# g_all = gridExtra::grid.arrange(g_resid_hist,g_qq,g_scatter_smooth,g_zero_hist,nrow=2)
g_all = gridExtra::grid.arrange(g_resid_hist,g_qq,g_zero_hist,
                                g_scatter_smooth,g_scatter_smooth3,g_scatter_smooth2,
                                ncol=2)

# ggsave(g_all, file="bestmodel_diagnosis.png",unit="mm",width=240,height=150,dpi=600)
# ggsave(g_all, file="bestmodel_diagnosis.png",unit="mm",width=240,height=150,dpi=600)
ggsave(g_all, file="bestmodel_diagnosis_rev.png",unit="mm",width=240,height=240,dpi=600)

(stats::cov2cor(bestmodel$rep$cov.fixed))

(time_peak = bestmodel$rep_summary[rownames(bestmodel$rep_summary) == "b_par",1]*bestmodel$time_sd+bestmodel$time_mean)
(time_peak_sd = bestmodel$rep_summary[rownames(bestmodel$rep_summary) == "b_par",2]*bestmodel$time_sd)


## index parameter plot ----

bestmodel$rep_summary %>% rownames %>% unique

par_pred = bestmodel$index %>%
  mutate(a_se = bestmodel$rep_summary[rownames(bestmodel$rep_summary)=="a_par",2],
         b_se = bestmodel$rep_summary[rownames(bestmodel$rep_summary)=="b_par",2],
         c_se = bestmodel$rep_summary[rownames(bestmodel$rep_summary)=="c_par",2]) %>%
  mutate(a = a_par-2*log(bestmodel$time_sd),
         b = b_par*bestmodel$time_sd+bestmodel$time_mean,
         b_se = b_se*bestmodel$time_sd) %>%
  rename(c = c_par) %>%
  dplyr::select(Year,a,b,c,a_se,b_se,c_se)
par_pred_mean = par_pred %>% dplyr::select(Year,a,b,c) %>%
  gather(key=parameter,value=value,-Year)
par_pred_se = par_pred %>% dplyr::select(Year,a_se,b_se,c_se) %>%
  rename(a=a_se,b=b_se,c=c_se) %>%
  gather(key=parameter,value=SE,-Year)
CI = 0.8
par_pred2 = full_join(par_pred_mean,par_pred_se) %>%
  mutate(upper = qnorm((1+CI)/2,value,SE),lower=qnorm((1-CI)/2,value,SE))

par_pred2_a = par_pred2 %>% dplyr::filter(parameter=="a")
par_pred2_b = par_pred2 %>% dplyr::filter(parameter=="b")
par_pred2_c = par_pred2 %>% dplyr::filter(parameter=="c")

g_a = ggplot(data=par_pred2_a,aes(x=Year)) +
  geom_ribbon(aes(ymax=upper,ymin=lower),alpha=0.4)+
  geom_path(aes(y=value),size =1)+
  geom_point(aes(y=value),size=2)+
  theme_bw(base_size=11)+ylab("a")+
  # geom_hline(yintercept=bestmodel$opt$par["a_mean"]-2*log(bestmodel$time_sd),linetype="dashed",colour="blue",
  # size=1)+
  ggtitle("(a) Pulse width")
g_a

g_b = ggplot(data=par_pred2_b,aes(x=Year)) +
  geom_ribbon(aes(ymax=upper,ymin=lower),alpha=0.4)+
  geom_path(aes(y=value),size =1)+
  geom_point(aes(y=value),size=2)+
  theme_bw(base_size=11)+ylab("b")+
  ggtitle("(b) Peak timing")+
  # geom_hline(yintercept=bestmodel$opt$par["b_mean"]*bestmodel$time_sd+bestmodel$time_mean,linetype="dashed",colour="blue",
  #            size=1)+
  ylim(min(par_pred2_b$lower)-3,max(par_pred2_b$upper)+3)
g_b

g_c = ggplot(data=par_pred2_c,aes(x=Year)) +
  geom_ribbon(aes(ymax=upper,ymin=lower),alpha=0.4)+
  geom_path(aes(y=value),size=1)+
  geom_point(aes(y=value),size=2)+
  theme_bw(base_size=11)+ylab("c")+
  ggtitle("(c) Peak size")
g_c

g_abc = gridExtra::grid.arrange(g_a,g_b,g_c,nrow=1)
ggsave(g_abc, file="parameter_change.png",dpi=600,width=240,height=80,unit="mm")


index = bestmodel$index %>%
  dplyr::select(Year,Estimate,CV) %>%
  rename("Standardized index" = Estimate)
index_scaled = bestmodel$index_scaled %>%
  dplyr::select(Year,Estimate,CV) %>%
  rename("Standardized index (scaled)" = Estimate, "CV (scaled)" = CV)

index = full_join(index,index_scaled,by="Year")

nominal = bestmodel$input$dat %>% group_by(Year) %>% 
  summarise(Max = max(Catch/Effort),Mean = mean(Catch/Effort))

index = full_join(nominal,index,by="Year")
index = index[,c(1:4,6,5,7)]

index2 = data.frame(Year = index$Year,
                    Nominal_max = str_c(round(index$Max,2)," (",round(index$Max/mean(index$Max),2),")",sep=""),
                    Nominal_mean = str_c(round(index$Mean,2)," (",round(index$Mean/mean(index$Mean),2),")",sep=""),
                    Standardized_index = str_c(round(index$`Standardized index`,2)," (",round(index$`Standardized index (scaled)`,2),")",sep=""),
                    CV = str_c(round(index$CV,2)," (",round(index$`CV (scaled)`,2),")",sep=""))

write.csv(index2, file="index_table.csv")

g_trend = plot_spm_trend2(bestmodel,CI=CI,legend_position=c(0.6,0.9))
g_trend$graph
ggsave(g_trend$graph, file="year_trend.png",dpi=600,width=180,height=100,unit="mm")

# bestmodel$opt$par["log_sigma_c"] %>% exp()

### retrospective analysis of standardized index ----

load("bestmodel.rda")

nretro = 10
# i=1
retro_index_dat = bestmodel$index %>% mutate(retro_year = 0)
retro_index_scaled_dat = bestmodel$index_scaled %>% mutate(retro_year = 0)
retro_list <- list()
for (i in 1:nretro) {
  if (i==1) {
    input = bestmodel$input
    # input$p0_list <- bestmodel$par_list
  }else{
    # input$p0_list <- model$par_list
  } 
  input$dat <- input$dat %>% dplyr::filter(Year < max(Year))
  retromodel = do.call(spm,input)
  retro_list[[i]] <- retromodel
  retro_index_dat = retro_index_dat %>% 
    bind_rows(retromodel$index %>% mutate(retro_year = i))
  retro_index_scaled_dat = retro_index_scaled_dat %>% 
    bind_rows(retromodel$index %>% mutate(retro_year = i))
}

retro_list[[10]]$rep

# colnames(retro_index_dat)

retro_index_dat_wider = retro_index_dat %>% 
  dplyr::select(Year,Estimate,retro_year) %>% 
  pivot_wider(names_from="retro_year",values_from="Estimate")

write.csv(retro_index_dat,file="retro_index_dat.csv",row.names=FALSE)
write.csv(retro_index_dat_wider,file="retro_index_dat_wider.csv",row.names=FALSE)

retro_index_dat = retro_index_dat %>%
  mutate(retro_id = as.character(retro_year))

rho_index = sapply(1:nretro,function(i) {
  tmp_dat = retro_index_dat_wider %>% select(Year,as.character(0),as.character((i))) %>% na.omit()
  as.numeric(tmp_dat[nrow(tmp_dat),c(3)]/tmp_dat[nrow(tmp_dat),c(2)]-1)
})
rho_index_dat = data.frame(Estimate=rho_index %>% mean(),Year=bestmodel$input$dat$Year %>% max()) %>%
  mutate(label=sprintf("rho == %.2f",Estimate)) %>%
  mutate(y = retro_index_dat$Estimate %>% max)

g_index_retro = ggplot(data=retro_index_dat)+
  geom_path(size=1,
            aes(x=Year,y=Estimate,colour=retro_id,linetype=retro_id))+
  scale_color_discrete(name="")+
  scale_linetype_discrete(name="")+
  theme_bw(base_size=12)+theme(legend.position="none")+
  geom_text(data=rho_index_dat,parse=TRUE,aes(x=Year,y=y,label=label,hjust=1,vjust=1))+
  ylab("Standardized index value")
g_index_retro

ggsave(g_index_retro,filename="index_retro.png",dpi=600,unit="mm",
       height=100,width=160)


### 仮に２０２１年のデータがすべてゼロだった場合

Dat3 = Dat2 %>% mutate(Catch = ifelse(Year==2021,0,Catch))

input_tmp = bestmodel$input
input_tmp$dat <- Dat3
bestmodel_zero = do.call(spm,input_tmp)
bestmodel_zero$index
