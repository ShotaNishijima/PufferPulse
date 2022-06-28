
#####################################################################
##### Analyzing VPA with short-term forecasting and hindcasting #####
#####################################################################

# Set working directory and read libraries ----

setwd("~/git/PufferPulse")

library(tidyverse)
# devtools::install_github("ichimomo/frasyr@dev")
library(frasyr)
library(mgcv)

# Read source code and best model ----

source("surfnet_peak_modeling.R")

load("bestmodel.rda")
dat = bestmodel$input$dat

nominal = dat %>% group_by(Year) %>%
  summarise(Max = max(Catch/Effort),Mean = mean(Catch/Effort))

# VPA analysis ----

#VPA data
caa <- read.csv("data/caa.csv",row.names=1)
waa <- read.csv("data/waa.csv",row.names=1)
maa <- read.csv("data/maa.csv",row.names=1)
M <- read.csv("data/M.csv",row.names=1)
cpue0 <- read.csv("data/cpue_age1.csv",row.names=1)

dat0 <- data.handler(caa, waa, maa, cpue0, M)

tf.mat <- matrix(NA,ncol=ncol(caa),nrow=nrow(caa))
tf.mat[1,(ncol(caa)-2):(ncol(caa)-1)] <- 1 

tf_year = as.numeric(rev(colnames(dat0$caa))[1])

## without age-0 index ----
vout <- vpa(dat0,
            tf.year=(tf_year-3):(tf_year-1),
            Pope=TRUE,
            fc.year=(tf_year-2):(tf_year-0),
            alpha=1,
            p.init=0.5,
            tune=TRUE,
            sel.update=TRUE,
            abund=c("N"),
            min.age=c(1),
            max.age=c(1),
            term.F = "max",
            tf.mat=tf.mat) 

vout$faa[,as.character((tf_year-4):(tf_year))]
vout$naa[,as.character((tf_year-4):(tf_year))]
vout$faa[,as.character(tf_year)]
dat0$caa[,as.character(tf_year)]


## Processing age-0 index data ----

# index_age0 = read.csv("../index_bestmodel.csv",header=T)
index_age0 = bestmodel$index

# index_age0 = read.csv("../index_bestmodel.csv",header=T)
rec_est0 = read.csv(paste0("data/Recruit_estimate",tf_year,".csv"), header=T) 
rec_est = rec_est0 %>% dplyr::filter(Year %in% unique(dat$Year)) %>%
  dplyr::mutate(Prop_nat = Natural/Total)

index_age0 = full_join(index_age0,rec_est) %>% 
  full_join(nominal) %>% 
  dplyr::select(Year,Estimate,Prop_nat,Mean,Max) %>%
  # dplyr::mutate(Mean2 = Mean+0.5*min(Mean[Mean>0]),
  #               Max2 = Max+0.5*min(Max[Max>0]))
  dplyr::mutate(Mean2 = Mean+1*min(Mean[Mean>0]),
                Max2 = Max+1*min(Max[Max>0]))

index_age0 = index_age0 %>%
  tidyr::replace_na(list(Prop_nat=1)) %>%
  dplyr::mutate(Standardized=Estimate,
                Nominal_max=Max2,Nominal_mean=Mean2)

dat1 <- dat2 <- dat3 <- dat0

# Standardized
dat1$index[2,colnames(dat1$index) %in% as.character(index_age0$Year)] <- index_age0 %>% filter(Year %in% as.numeric(colnames(dat0$index))) %>%
  dplyr::select(Standardized) %>% unlist() %>% as.numeric()

# Nominal Max 
dat2$index[2,colnames(dat2$index) %in% as.character(index_age0$Year)] <- index_age0 %>% filter(Year %in% as.numeric(colnames(dat0$index))) %>%
  dplyr::select(Nominal_max) %>% unlist() %>% as.numeric()

# Nominal Mean 
dat3$index[2,colnames(dat3$index) %in% as.character(index_age0$Year)] <- index_age0 %>% filter(Year %in% as.numeric(colnames(dat0$index))) %>%
  dplyr::select(Nominal_mean) %>% unlist() %>% as.numeric()

row.names(dat1$index) <- row.names(dat2$index) <- row.names(dat3$index) <- c("age1","age0")

remove.abund = dat1$index
remove.abund[1,] <- 0
remove.abund[2,] <- rec_est0$Additive*exp(0.125) #Transition from October -> April


## vpa analysis with age-0 index ----

#with standardized
vout1 <- vpa(dat1,
             tf.year=(tf_year-3):(tf_year-1),
             Pope=TRUE,
             fc.year=(tf_year-2):tf_year,
             alpha=1,
             p.init=c(0.8),
             tune=TRUE,
             sel.update=TRUE,
             abund=c("N","N"),
             min.age=c(1,0),
             max.age=c(1,0),
             term.F = "max",
             add.p.est = 0,
             add.p.ini=0.1,
             est.method="ls", #change to LS 2022/06/15
             remove.abund=remove.abund) 

vout$input$fc.year
dat1$index
vout1$term.f
vout1$faa[,as.character((tf_year-4):tf_year)]
vout1$sigma
vout1$pred.index

# nominal max
input_tmp = vout1$input
input_tmp$dat = dat2
input_tmp$add.p.ini=0.2
vout2 = do.call(vpa, input_tmp)
vout2$sigma
vout2$term.f
vout2$faa[,as.character((tf_year-4):tf_year)]

# nominal mean
input_tmp = vout1$input
input_tmp$dat = dat3
input_tmp$add.p.ini=0.2
vout3 = do.call(vpa, input_tmp)
vout3$sigma
vout3$term.f
vout3$faa[,as.character((tf_year-4):tf_year)]

cbind(vout1$sigma,vout2$sigma,vout3$sigma)

#the best fit with standardized

## residual plot ----

residual1 = plot_residual_vpa(vout1,index_name=c("Age 1", "Age 0"),plot_year=1995:2020)
# ggsave(residual1$year_resid + ylab("Residual"),filename="res/year_resid.png",dpi=600,unit="mm",height=80,width=180)
# ggsave(residual1$fitting_Index,filename="res/fitting_index.png",dpi=600,unit="mm",height=80,width=180)
# ggsave(residual1$abund_Index + xlab("Abundance"),filename="res/abund_index.png",dpi=600,unit="mm",height=80,width=180)


resid1_age0 = residual1$gg_data %>% filter(Index_Label=="Age 0") %>%
  dplyr::select(year,resid) %>% na.omit()


m1 = gam(resid~s(year),data=resid1_age0)
summary(m1)

g1 = ggplot(data=resid1_age0,aes(year,resid))+
  geom_point(size=2)+theme_bw()+
  geom_hline(yintercept=0)+
  geom_smooth(method="gam")

g1



residual2 = plot_residual_vpa(vout2,index_name=c("Age 1", "Age 0"))
residual2$year_resid
residual2$fitting_Index

resid2_age0 = residual2$gg_data %>% filter(Index_Label=="Age 0") %>%
  dplyr::select(year,resid) %>% na.omit()

m2 = gam(resid~s(year),data=resid2_age0)
summary(m2)

g2 = ggplot(data=resid2_age0,aes(year,resid))+
  geom_point(size=2)+theme_bw()+
  geom_hline(yintercept=0)+
  geom_smooth(method="gam")

g2

residual3 = plot_residual_vpa(vout3,index_name=c("Age 1", "Age 0"))
residual3$year_resid
residual3$fitting_Index

resid3_age0 = residual3$gg_data %>% filter(Index_Label=="Age 0") %>%
  dplyr::select(year,resid) %>% na.omit()

m3 = gam(resid~s(year),data=resid3_age0)
summary(m3)

g3 = ggplot(data=resid3_age0,aes(year,resid))+
  geom_point(size=2)+theme_bw()+
  geom_hline(yintercept=0)+
  geom_smooth(method="gam")

g3

resid1_age0 = residual1$gg_data %>% filter(Index_Label=="Age 0") %>%
  dplyr::mutate(Index = "Standardized")

resid2_age0 = residual2$gg_data %>% filter(Index_Label=="Age 0") %>%
  dplyr::mutate(Index = 'Nominal mean')

resid3_age0 = residual3$gg_data %>% filter(Index_Label=="Age 0") %>%
  dplyr::mutate(Index = 'Nominal max')

resid_age0 = full_join(resid1_age0,resid2_age0) %>% full_join(resid3_age0) %>%
  na.omit() %>%
  rename(Year=year)

g_fit = ggplot(resid_age0,aes(x=Year)) +
  geom_path(aes(y=pred),size=1.2) +
  geom_point(aes(y=obs),size=2) +
  ylab("Index")+
  facet_wrap(vars(Index),ncol=3,scales="free_y")+
  theme_bw(base_size=16)+
  ggtitle("(a) Prediction vs. Observation")
g_fit


resid_label = resid_age0 %>% 
  mutate(y = min(resid)) %>%
  dplyr::filter(Year==min(Year)) %>%
  dplyr::select(Index,Year,y,sigma) %>%
  mutate(label = sprintf("sigma == %.2f",sigma))

p_label = resid_label %>%
  mutate(y = min(resid_age0$resid)) %>%
  mutate(p = c(summary(m1)[["s.table"]][1,4],summary(m2)[["s.table"]][1,4],summary(m3)[["s.table"]][1,4])) %>%
  mutate(label= sprintf("P == %.2f",p)) %>%
  mutate(y = 0.9*y)

g_resid = ggplot(resid_age0,aes(x=Year)) +
  geom_point(aes(y=resid),size=2) +
  ylab("Residual")+
  facet_wrap(vars(Index),ncol=3,scales="fixed")+
  theme_bw(base_size=16)+
  geom_smooth(aes(y=resid),method="gam",size=1.2)+
  geom_hline(yintercept=0)+
  geom_text(data=resid_label,parse=TRUE,aes(y=y,label=label,hjust=0,vjust=0.5),size=5)+
  geom_text(data=p_label,parse=TRUE,aes(y=y,label=label,hjust=0,vjust=0),size=5)+
  ggtitle("(b) Residual trend")
g_resid


g_both = gridExtra::grid.arrange(g_fit,g_resid,nrow=2)
g_both

ggsave(g_both,file="res/VPA_fit.png",dpi=600,unit="mm",height=180,width=240)


## retrospective analysis ----

vout_retro = do_retrospective_vpa(vout,n_retro=10,plot_year = (tf_year-14):tf_year)
vout_retro$mohn_rho
vout_retro$graph

grid.add.ini <- seq(0.1,1.5,by=0.05)
vout1$input$p.init <- 1.5
vout1_retro = do_retrospective_vpa(vout1,n_retro = 10,grid_add_ini = grid.add.ini,plot_year = (tf_year-14):tf_year)
vout1_retro$mohn_rho
vout1_retro$graph

# grid.add.ini <- seq(0.1,1.5,by=0.05)
vout2$input$p.init <- 1.5
vout2_retro = do_retrospective_vpa(vout2,n_retro=10,grid_add_ini = grid.add.ini,plot_year=(tf_year-14):tf_year)
vout2_retro$mohn_rho
vout2_retro$graph

vout3$input$p.init <- 1.5
grid.add.ini <- seq(0.1,1.5,by=0.05)
# grid.init = seq(0.1,2,by=0.1)
vout3_retro = do_retrospective_vpa(vout3,n_retro=10,grid_init = NULL,
                                   grid_add_ini = grid.add.ini,plot_year=(tf_year-14):tf_year)
vout3_retro$mohn_rho
vout3_retro$result$retro.f


## Compare estimates ---- 

vpa_list = list(vout,vout1,vout2,vout3)

g_vpa = plot_vpa(vpa_list,what.plot=c("Recruitment","biomass","SSB","U"),
                 vpaname=c("w/o age0-index","Standardized","Nominal mean","Nominal max"))
labeli = as_labeller(c('Recruitment'="Recruitment",'biomass'="Biomass",'SSB'="SSB",'U'="U"))

vpa_dat=g_vpa$data %>%
  dplyr::filter(year>2014) %>%
  dplyr::mutate(value = if_else(value>1,value/1000,value)) %>%
  dplyr::mutate(id=factor(id,levels=c("Standardized","Nominal max","Nominal mean","w/o age0-index")))

g_vpa2 = ggplot(data=vpa_dat,aes(x=year,y=value,group=id,linetype=id,colour=id,shape=id))+
  geom_path(size=1)+geom_point(size=1.5)+
  facet_wrap(vars(stat),scales="free_y",labeller = labeli)+
  scale_colour_brewer(palette="Set1",name="")+
  scale_shape(solid=TRUE,name="")+
  scale_linetype(name="")+
  theme_bw()+ylab("")+ylim(0,NA)+xlab("Year")+
  theme_SH(legend.position="bottom")
g_vpa2

ggsave(g_vpa2,filename="res/compare_vpa.png",unit="mm",dpi=600,height=120,width=180)


## retrospective forecasting ----

rec_additive = rec_est0$Additive*exp(0.125)
# vpares =vout
# index_value = NULL

# without age0 index
vout_forecast = one_year_forecast(vout, index_value=NULL,rec_additive=rec_additive)
vout_forecast = vout_forecast %>% mutate(retro_year = 0)

max_year = max(as.numeric(colnames(vout$naa)))
for(i in 1:length(vout_retro$result$Res)) {
  rec_additive2 = pull(filter(rec_est0,Year<=max_year-i),Additive)
  df = one_year_forecast(vout_retro$result$Res[[i]],index_value=NULL,rec_additive=rec_additive2) %>% dplyr::mutate(retro_year=i)
  vout_forecast = full_join(vout_forecast, df)
}

# with standardized index
vout1_forecast = one_year_forecast(vout1, index_value=rev(index_age0$Estimate)[1],rec_additive=rec_additive)
vout1_forecast = vout1_forecast %>% mutate(retro_year = 0)

for(i in 1:length(vout1_retro$result$Res)) {
  rec_additive2 = pull(filter(rec_est0,Year<=max_year-i),Additive)
  df = one_year_forecast(vout1_retro$result$Res[[i]],index_value=rev(index_age0$Estimate)[1+i],rec_additive=rec_additive2) %>% dplyr::mutate(retro_year=i)
  vout1_forecast = full_join(vout1_forecast, df)
}


# with nominal max
vout2_forecast = one_year_forecast(vout2, index_value=rev(index_age0$Max2)[1],rec_additive=rec_additive)
vout2_forecast = vout2_forecast %>% mutate(retro_year = 0)

for(i in 1:length(vout2_retro$result$Res)) {
  rec_additive2 = pull(filter(rec_est0,Year<=max_year-i),Additive)
  df = one_year_forecast(vout2_retro$result$Res[[i]],index_value=rev(index_age0$Max2)[1+i],rec_additive=rec_additive2) %>% dplyr::mutate(retro_year=i)
  vout2_forecast = full_join(vout2_forecast, df)
}

# with nominal mean
vout3_forecast = one_year_forecast(vout3, index_value=rev(index_age0$Mean2)[1],rec_additive=rec_additive)
vout3_forecast = vout3_forecast %>% mutate(retro_year = 0)

for(i in 1:length(vout3_retro$result$Res)) {
  rec_additive2 = pull(filter(rec_est0,Year<=max_year-i),Additive)
  df = one_year_forecast(vout3_retro$result$Res[[i]],index_value=rev(index_age0$Mean2)[1+i],rec_additive=rec_additive2) %>% dplyr::mutate(retro_year=i)
  vout3_forecast = full_join(vout3_forecast, df)
}

vout_mohn = get_mohn_graph(vout_forecast)
g0 = vout_mohn$gg + ggtitle("Without age0-index")
g0
ggsave(g0,file="res/retro_without_age0.png",unit="mm",dpi=600,height=150,width=240)

vout1_mohn = get_mohn_graph(vout1_forecast)
g1 = vout1_mohn$gg + ggtitle("With standardized index")
g1
ggsave(g1,file="res/retro_standardized.png",unit="mm",dpi=600,height=150,width=240)

vout2_mohn = get_mohn_graph(vout2_forecast)
g2 = vout2_mohn$gg + ggtitle("With nominal-max index")
g2
ggsave(g2,file="res/retro_nominal-max.png",unit="mm",dpi=600,height=150,width=240)

vout3_mohn = get_mohn_graph(vout3_forecast)
g3 = vout3_mohn$gg + ggtitle("With nominal-mean index")
g3
ggsave(g3,file="res/retro_nominal-mean.png",unit="mm",dpi=600,height=150,width=240)

rf_dat = vout_mohn$dat_longer %>% mutate(Model = "No age0-index") %>%
  bind_rows(vout1_mohn$dat_longer %>% mutate(Model = "Standardized")) %>%
  bind_rows(vout2_mohn$dat_longer %>% mutate(Model = "Nominal-max")) %>%
  bind_rows(vout3_mohn$dat_longer %>% mutate(Model = "Nominal-mean")) %>%
  filter(Variable != "F") %>%
  filter(Year > (max(.$Year)-16)) %>%
  mutate(Model_f = factor(Model,levels=c("No age0-index","Nominal-max","Nominal-mean","Standardized"))) %>%
  mutate(retro_id = as.character(retro_year))

rho_dat = rownames_to_column(data.frame(value=vout_mohn$rho),"Variable") %>% mutate(Model = "No age0-index") %>%
  bind_rows(rownames_to_column(data.frame(value=vout1_mohn$rho),"Variable") %>% mutate(Model = "Standardized")) %>%
  bind_rows(rownames_to_column(data.frame(value=vout2_mohn$rho),"Variable") %>% mutate(Model = "Nominal-max")) %>%
  bind_rows(rownames_to_column(data.frame(value=vout3_mohn$rho),"Variable") %>% mutate(Model = "Nominal-mean")) %>%
  filter(Variable != "F" & Variable != "Number") %>%
  mutate(Year = max(vout1_mohn$dat_longer$Year)) %>%
  full_join(rf_dat %>% group_by(Variable) %>%
              summarise(y = max(Value))) %>% 
  mutate(label=sprintf("rho == %.2f",value)) %>%
  mutate(Variable=factor(Variable,levels=c("Recruit","Biomass","SSB","U","Catch"))) %>%
  mutate(Model_f = factor(Model,levels=c("No age0-index","Nominal-max","Nominal-mean","Standardized")))

g_retro_all = ggplot(data=rf_dat)+
  geom_path(size=0.8,aes(x=Year,y=Value,colour=retro_id,linetype=retro_id))+
  scale_colour_discrete(name="")+
  theme_bw(base_size=10)+theme(legend.position="none")+
  geom_text(data=rho_dat,parse=TRUE,size=3,aes(x=Year,y=y,label=label,hjust=1,vjust=1))+
  facet_grid(rows=vars(Variable),cols=vars(Model_f),scales="free_y")+
  ylab("")+ylim(0,NA)
g_retro_all  

rf_dat_full = rf_dat %>% filter(retro_id=="0")
rf_dat_others = rf_dat %>% filter(retro_id!="0")

rf_dat_full_point = rf_dat_full %>% ungroup() %>% filter(Year == max(Year))

rf_dat_others_point = rf_dat_others %>% ungroup() %>% group_by(retro_id) %>% 
  filter(Year == max(Year))

g_retro_all2 = ggplot(data=rf_dat_others)+
  geom_path(size=0.8,aes(x=Year,y=Value,colour=retro_id))+
  geom_point(data=rf_dat_others_point,aes(x=Year,y=Value,colour=retro_id),size=1.5)+
  geom_path(data=rf_dat_full,aes(x=Year,y=Value),size=0.8,colour="black")+
  geom_point(data=rf_dat_full_point,aes(x=Year,y=Value),size=1.5,colour="black")+
  scale_colour_discrete(name="")+
  theme_bw(base_size=10)+theme(legend.position="none")+
  geom_text(data=rho_dat,parse=TRUE,size=3,aes(x=Year,y=y,label=label,hjust=1,vjust=1))+
  facet_grid(rows=vars(Variable),cols=vars(Model_f),scales="free_y")+
  ylab("")+ylim(0,NA)
g_retro_all2  

ggsave(g_retro_all,filename="res/retro_all_models.png",dpi=600,unit="mm",
       height=160,width=240)

ggsave(g_retro_all2,filename="res/retro_all_models_rev.png",dpi=600,unit="mm",
       height=160,width=240)

dat_forecast = vout_mohn$dat_plot %>% mutate(Index = "w/o age0-index") %>%
  full_join(vout1_mohn$dat_plot %>% mutate(Index = "Standardized")) %>%
  full_join(vout2_mohn$dat_plot %>% mutate(Index = "Nominal max")) %>%
  full_join(vout3_mohn$dat_plot %>% mutate(Index = "Nominal mean"))

dat_forecast = dat_forecast %>% filter(retro_year=='0') %>%
  mutate(Index = factor(Index,levels=c("Standardized","Nominal mean","Nominal max","w/o age0-index"))) %>%
  filter(Year>2014)

g_forecast = ggplot(data= dat_forecast, aes(x=Year,y=Value,linetype=Index,shape=Index,colour=Index))+
  geom_point(size=1.5)+
  geom_path(size=1)+
  facet_wrap(vars(Variable),ncol=3,scales="free_y")+
  theme_bw(base_size=16)+
  scale_colour_brewer(palette="Set1",name="")+
  scale_linetype(name="")+
  scale_shape(name="")+
  theme(legend.position="bottom",legend.text=element_text(size=16))+
  ylim(0,NA)+ylab("")
g_forecast

# filter(dat_forecast,Index=="Standardized" & Variable=="Recruit") %>% View()

ggsave(g_forecast,file="res/compare_forecast.png",dpi=600,unit="mm",height=180,width=240)


out.vpa(vout,filename="res/Vpa_No_Age0_Index")
out.vpa(vout1,filename="res/Vpa_Age0_Index_standardized")
out.vpa(vout2,filename="res/Vpa_Age0_Index_nominal-max")
out.vpa(vout3,filename="res/Vpa_Age0_Index_nominal-mean")

### output selectivity ---- 

dat_saa = vout1$faa %>% as.data.frame() %>% rownames_to_column(var="age") %>%
  pivot_longer(cols=-age,names_to = "year",values_to = "faa") %>%
  mutate(Age = factor(age,levels=rev(c("0","1","2","3"))),Year=as.numeric(year)) %>%
  dplyr::select(Age,Year,faa) %>%
  arrange(Year,Age) %>%
  group_by(Year) %>%
  mutate(sumF = sum(faa)) %>%
  mutate(saa = faa/sumF) %>% ungroup
  

(g_saa = ggplot(filter(dat_saa,Age!="3"),aes(x=Year,y=saa,fill=Age))+
  geom_bar(stat="identity",position="fill")+
    scale_fill_brewer(palette="Accent")+
    theme_bw(base_size=12)+
  ylim(0,1)+ylab("Selectivity at age")
)

ggsave(g_saa,filename="res/selectivity.png",dpi=600,unit="mm",height=100,width=150)


## hindcast cross validation ----

## one year forecast
index_age1 = vout$input$dat$index %>% t() %>% as.data.frame() %>%
  rownames_to_column(., var = "rowname") %>%
  mutate(Year = as.numeric(rowname)) %>% 
  na.omit() %>% dplyr::select(-rowname) %>%
  rename(Observation = age1)


# no age 0 index
Q = vout_retro$result$Res %>% map_dbl(., function(x) x$q[1])
vout_cv = vout_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & Year == max(Year) & retro_year> 0) %>%
  ungroup() %>%
  mutate(q = Q) %>%
  full_join(index_age1) %>%
  mutate(Prediction = q*naa)

Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout_cv_full = vout_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & retro_year> 0) %>%
  ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa)

vout_MASE = data.frame(denominator=mean(abs(diff(pull(vout_cv %>% na.omit(),Observation)))),
                       numerator = mean(abs(pull(vout_cv %>% na.omit() %>% mutate(Error = Observation-Prediction),Error)))) %>%
  mutate(MASE = numerator/denominator)

# standardized index
Q = vout1_retro$result$Res %>% map_dbl(., function(x) x$q[1])
vout1_cv = vout1_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & Year == max(Year) & retro_year> 0) %>%
  ungroup() %>%
  mutate(q = Q) %>%
  full_join(index_age1) %>%
  mutate(Prediction = q*naa)


Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout1_cv_full = vout1_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & retro_year> 0) %>%
  ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa)

vout1_MASE = data.frame(denominator=mean(abs(diff(pull(vout1_cv %>% na.omit(),Observation)))),
                       numerator = mean(abs(pull(vout1_cv %>% na.omit() %>% mutate(Error = Observation-Prediction),Error)))) %>%
  mutate(MASE = numerator/denominator)

# sd(vout_cv$q,na.rm=TRUE)/mean(vout_cv$q,na.rm=TRUE)
# sd(vout1_cv$q,na.rm=TRUE)/mean(vout1_cv$q,na.rm=TRUE)
# matplot(cbind(vout_cv$q,vout1_cv$q) %>% na.omit())

# nominal max index
Q = vout2_retro$result$Res %>% map_dbl(., function(x) x$q[1])
vout2_cv = vout2_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & Year == max(Year) & retro_year> 0) %>%
  ungroup() %>%
  mutate(q = Q) %>%
  full_join(index_age1) %>%
  mutate(Prediction = q*naa)

Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout2_cv_full = vout2_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & retro_year> 0) %>%
  ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa)

vout2_MASE = data.frame(denominator=mean(abs(diff(pull(vout2_cv %>% na.omit(),Observation)))),
                        numerator = mean(abs(pull(vout2_cv %>% na.omit() %>% mutate(Error = Observation-Prediction),Error)))) %>%
  mutate(MASE = numerator/denominator)


# nominal mean index
Q = vout3_retro$result$Res %>% map_dbl(., function(x) x$q[1])
vout3_cv = vout3_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & Year == max(Year) & retro_year> 0) %>%
  ungroup() %>%
  mutate(q = Q) %>%
  full_join(index_age1) %>%
  mutate(Prediction = q*naa)

Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout3_cv_full = vout3_forecast %>% group_by(retro_year) %>%
  filter(Age==1 & retro_year> 0) %>%
  ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa)


vout3_MASE = data.frame(denominator=mean(abs(diff(pull(vout3_cv %>% na.omit(),Observation)))),
                        numerator = mean(abs(pull(vout3_cv %>% na.omit() %>% mutate(Error = Observation-Prediction),Error)))) %>%
  mutate(MASE = numerator/denominator)

dat_MASE = vout_MASE %>% mutate(Index = "No age0-index") %>%
  full_join(vout1_MASE %>% mutate(Index = "Standardized")) %>%
  full_join(vout2_MASE %>% mutate(Index = "Nominal-max")) %>%
  full_join(vout3_MASE %>% mutate(Index = "Nominal-mean")) %>%
  mutate(Index = factor(Index,levels=c("No age0-index","Nominal-max","Nominal-mean","Standardized")))

dat_MASE = dat_MASE %>%
  mutate(label = sprintf("MASE == %.2f",MASE))

dat_hindCV = vout_cv_full %>% mutate(Index = "No age0-index") %>%
  full_join(vout1_cv_full %>% mutate(Index = "Standardized")) %>%
  full_join(vout2_cv_full %>% mutate(Index = "Nominal-max")) %>%
  full_join(vout3_cv_full %>% mutate(Index = "Nominal-mean")) %>%
  mutate(Index = factor(Index,levels=c("No age0-index","Nominal-max","Nominal-mean","Standardized"))) %>%
  na.omit() %>%
  mutate(retro_id = as.character(retro_year)) %>%
  filter(Year>2005) %>% 
  arrange(Index,Year)

dat_hindCV_full = dat_hindCV %>% filter(retro_year==1)
# dat_hindCV2_full %>% View
dat_hindCV_point = dat_hindCV %>% group_by(retro_id) %>%
  filter(Year==max(Year)) %>% ungroup()

(g_hindCV = ggplot(data=dat_hindCV,aes(x=Year))+
    geom_path(aes(y=Prediction,colour=retro_id),size=0.8)+
    geom_point(data=dat_hindCV_point,aes(y=Prediction,colour=retro_id),size=2)+
    geom_path(data=dat_hindCV_full,aes(y=Observation),colour="black",size=0.8)+
    scale_colour_discrete(name="")+
    geom_text(data=dat_MASE,parse=TRUE,size=4,aes(x=2020,y=220500,label=label,hjust=1,vjust=1))+
    facet_wrap(~Index)+
    theme_bw(base_size=12)+theme(legend.position="none")+
    scale_y_continuous(labels = scales::label_comma())+
    ylab("Age-1 index value")
)

ggsave(g_hindCV,file="res/hindcast_CV.png",dpi=600,unit="mm",height=120,width=180)


## two year forecast ----

Q = vout_retro$result$Res %>% map_dbl(., function(x) x$q[1])
Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout_cv_full2 = vout_forecast %>% group_by(retro_year) %>%
  filter(Year==max(Year) & Age==0) %>%
  mutate(naa = naa*exp(-0.25-faa),
         Age=Age+1,Year=Year+1)  %>% ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa) %>%
  full_join(.,vout_cv_full)

temp = vout_cv_full2 %>% 
  filter(retro_year>1) %>% 
  group_by(retro_year) %>% 
  filter(Year==max(Year)) %>% ungroup() %>%
  mutate(Error = Observation-Prediction)

vout_MASE2 = data.frame(numerator=mean(abs(temp$Error)),
                        denomator=mean(abs(diff(pull(index_age1 %>% filter(Year >= min(temp$Year)-4),Observation),lag=2)))) %>%
  mutate(MASE = numerator/denomator)

Q = vout1_retro$result$Res %>% map_dbl(., function(x) x$q[1])
Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout1_cv_full2 = vout1_forecast %>% group_by(retro_year) %>%
  filter(Year==max(Year) & Age==0) %>%
  mutate(naa = naa*exp(-0.25-faa),
         Age=Age+1,Year=Year+1)  %>% ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa) %>%
  full_join(.,vout1_cv_full)

temp = vout1_cv_full2 %>% 
  filter(retro_year>1) %>% 
  group_by(retro_year) %>% 
  filter(Year==max(Year)) %>% ungroup() %>%
  mutate(Error = Observation-Prediction)

vout1_MASE2 = data.frame(numerator=mean(abs(temp$Error)),
                         denomator=mean(abs(diff(pull(index_age1 %>% filter(Year >= min(temp$Year)-4),Observation),lag=2)))) %>%
  mutate(MASE = numerator/denomator)

Q = vout2_retro$result$Res %>% map_dbl(., function(x) x$q[1])
Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout2_cv_full2 = vout2_forecast %>% group_by(retro_year) %>%
  filter(Year==max(Year) & Age==0) %>%
  mutate(naa = naa*exp(-0.25-faa),
         Age=Age+1,Year=Year+1)  %>% ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa) %>%
  full_join(.,vout2_cv_full)

temp = vout2_cv_full2 %>% 
  filter(retro_year>1) %>% 
  group_by(retro_year) %>% 
  filter(Year==max(Year)) %>% ungroup() %>%
  mutate(Error = Observation-Prediction)

vout2_MASE2 = data.frame(numerator=mean(abs(temp$Error)),
                         denomator=mean(abs(diff(pull(index_age1 %>% filter(Year >= min(temp$Year)-4),Observation),lag=2)))) %>%
  mutate(MASE = numerator/denomator)

Q = vout3_retro$result$Res %>% map_dbl(., function(x) x$q[1])
Q_dat = data.frame(q = Q) %>% mutate(retro_year = 1:n())
vout3_cv_full2 = vout3_forecast %>% group_by(retro_year) %>%
  filter(Year==max(Year) & Age==0) %>%
  mutate(naa = naa*exp(-0.25-faa),
         Age=Age+1,Year=Year+1)  %>% ungroup() %>%
  full_join(index_age1) %>%
  full_join(Q_dat) %>%
  mutate(Prediction = q*naa) %>%
  full_join(.,vout3_cv_full)

temp = vout3_cv_full2 %>% 
  filter(retro_year>1) %>% 
  group_by(retro_year) %>% 
  filter(Year==max(Year)) %>% ungroup() %>%
  mutate(Error = Observation-Prediction)

vout3_MASE2 = data.frame(numerator=mean(abs(temp$Error)),
                         denomator=mean(abs(diff(pull(index_age1 %>% filter(Year >= min(temp$Year)-4),Observation),lag=2)))) %>%
  mutate(MASE = numerator/denomator)

dat_MASE2 = vout_MASE2 %>% mutate(Index = "No age0-index") %>%
  full_join(vout1_MASE2 %>% mutate(Index = "Standardized")) %>%
  full_join(vout2_MASE2 %>% mutate(Index = "Nominal-max")) %>%
  full_join(vout3_MASE2 %>% mutate(Index = "Nominal-mean")) %>%
  mutate(Index = factor(Index,levels=c("No age0-index","Nominal-max","Nominal-mean","Standardized")))

dat_MASE2 = dat_MASE2 %>%
  mutate(label = sprintf("MASE == %.2f",MASE))

dat_hindCV2 = vout_cv_full2 %>% mutate(Index = "No age0-index") %>%
  full_join(vout1_cv_full2 %>% mutate(Index = "Standardized")) %>%
  full_join(vout2_cv_full2 %>% mutate(Index = "Nominal-max")) %>%
  full_join(vout3_cv_full2 %>% mutate(Index = "Nominal-mean")) %>%
  mutate(Index = factor(Index,levels=c("No age0-index","Nominal-max","Nominal-mean","Standardized"))) %>%
  na.omit() %>%
  mutate(retro_id = as.character(retro_year)) %>%
  filter(Year>2005 & retro_year>1) %>% 
  arrange(Index,Year)

dat_hindCV2_full = dat_hindCV2 %>% filter(retro_year==2)
# dat_hindCV2_full %>% View
dat_hindCV2_point = dat_hindCV2 %>% group_by(retro_id) %>%
  filter(Year==max(Year)) %>% ungroup()

(g_hindCV2 = ggplot(data=dat_hindCV2,aes(x=Year))+
    geom_path(aes(y=Prediction,colour=retro_id),size=0.8)+
    geom_point(data=dat_hindCV2_point,aes(y=Prediction,colour=retro_id),size=2)+
    geom_path(data=dat_hindCV2_full,aes(y=Observation),colour="black",size=0.8)+
    scale_colour_discrete(name="")+
    geom_text(data=dat_MASE2,parse=TRUE,size=4,aes(x=2020,y=220500,label=label,hjust=1,vjust=1))+
    facet_wrap(~Index)+
    theme_bw(base_size=12)+theme(legend.position="none")+
    scale_y_continuous(labels = scales::label_comma())+
    ylab("Age-1 index value")
)

ggsave(g_hindCV2,file="res/hindcast_CV2.png",dpi=600,unit="mm",height=120,width=180)


### Figure h=1 and h=2 ----

(g_hindCV_both1 = ggplot(data=dat_hindCV,aes(x=Year))+
    geom_path(aes(y=Prediction,colour=retro_id),size=0.8)+
    geom_point(data=dat_hindCV_point,aes(y=Prediction,colour=retro_id),size=1.5)+
    geom_path(data=dat_hindCV_full,aes(y=Observation),colour="black",size=0.8)+
    scale_colour_discrete(name="")+
    geom_text(data=dat_MASE,parse=TRUE,size=3,aes(x=2020,y=240000,label=label,hjust=1,vjust=1))+
    facet_wrap(~Index,nrow=1)+
    theme_bw(base_size=10)+theme(legend.position="none")+
    scale_y_continuous(labels = scales::label_comma())+
    ylab("Age-1 index value")+
    ggtitle("(a) one-year ahead")
)

(g_hindCV_both2 = ggplot(data=dat_hindCV2,aes(x=Year))+
    geom_path(aes(y=Prediction,colour=retro_id),size=0.8)+
    geom_point(data=dat_hindCV2_point,aes(y=Prediction,colour=retro_id),size=1.5)+
    geom_path(data=dat_hindCV2_full,aes(y=Observation),colour="black",size=0.8)+
    scale_colour_discrete(name="")+
    geom_text(data=dat_MASE2,parse=TRUE,size=3,aes(x=2020,y=240000,label=label,hjust=1,vjust=1))+
    facet_wrap(~Index,nrow=1)+
    theme_bw(base_size=10)+theme(legend.position="none")+
    scale_y_continuous(labels = scales::label_comma())+
    ylab("Age-1 index value")+
    ggtitle("(b) two-year ahead")
  
)

g_hindCV_both =gridExtra::grid.arrange(g_hindCV_both1,g_hindCV_both2)

ggsave(g_hindCV_both,file="res/hindcast_CV_both.png",dpi=600,unit="mm",height=150,width=240)


## save forecasting results ----

dat_forecast_all = vout_forecast %>% mutate(Model = "No-age0-index") %>%
  bind_rows(vout1_forecast %>% mutate(Model = "Standardized") ) %>%
  bind_rows(vout2_forecast %>% mutate(Model = "Nominal-max") ) %>%
  bind_rows(vout3_forecast %>% mutate(Model = "Nominal-mean") )

write.csv(dat_forecast_all,file="res/dat_forecast_all.csv",row.names=FALSE)

### Fin ----