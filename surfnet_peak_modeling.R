


spm = function(dat,
               family="poisson", # or poisson or nbinom2
               a_key = 1, ## 0: fixed effect, 1: constant, 2: white noise, 3: AR1, 4: Random walk
               b_key = 1,
               c_key = 2,
               scaling_time = TRUE,
               # dat_tuningCPUE = NULL,
               p0_list = NULL,
               silent=FALSE,
               bias_correct=TRUE,
               w = NULL,
               use_MVN = 0,
               map_psi = 0:2,
               map_add = NULL,
               sigma_init = NULL) {
  
  argname <- ls()
  arglist <- lapply(argname,function(xx) eval(parse(text=xx)))
  names(arglist) <- argname
  
  Res = list()
  Res$input <- arglist
  
  Catch = dat$Catch
  Effort = dat$Effort
  iYear = dat$Year - min(dat$Year)
  if (scaling_time) {
    time_mean = dat$Time %>% mean()
    time_sd = dat$Time %>% sd()
    Res$time_mean = time_mean
    Res$time_sd = time_sd
    Time = (dat$Time - time_mean)/time_sd
    Time_grid = (min(dat$Time):max(dat$Time) - time_mean)/time_sd
  } else {
    Time = dat$Time
    Time_grid = min(dat$Time):max(dat$Time)
  }
  
  # family = "nbinom1"
  # family = "poisson"
  # family= "nbinom2"
  if (family == "poisson") family_tmb = 0
  if (family == "nbinom1") family_tmb = 1
  if (family == "nbinom2") family_tmb = 2
  
  if (is.null(w)) w = rep(1, length(Catch))
  
  data = list(Catch=Catch,Effort=Effort,iYear=iYear,Time=Time,
              a_key=a_key,b_key=b_key,c_key=c_key,family=family_tmb,
              Time_grid=Time_grid, w=w, Time_sd=time_sd,
              use_MVN=use_MVN
              # ,tuningCPUE_key=tuningCPUE_key,tuningCPUE=tuningCPUE,iYear_tuningCPUE=iYear_tuningCPUE
              )
  ny = max(iYear)+1
  
  if (is.null(p0_list)) {
    # glmmres = glmmTMB::glmmTMB(Catch ~ 1 + (1|as.factor(dat$Year)) + Time + I(Time^2),
    #                            offset = log(Effort), family = family)
    glmdata = data.frame(Catch=Catch,Effort=Effort,iYear = as.factor(iYear),Time=Time)
    # glmdata$iYear
    glmmres = glmmTMB::glmmTMB(Catch ~ 1 + (1|iYear) + Time + I(Time^2),
                               offset = log(Effort), family = family,data=glmdata)
    # summary(glmmres)
    coefs = coef(glmmres)$cond$'iYear'
    # coefs = coef(glmmres)$cond$'as.factor(dat$Year)'
    a_mean = log(-coefs[1,3]) 
    b_mean = coefs[1,2]
    c_mean = mean(coefs[,1])
    
    a_y = rep(a_mean,ny)
    b_y = rep(b_mean,ny)
    # c_y = rep(c_mean,ny)
    c_y = coefs[,1]
    
    log_sigma_a <- log_sigma_b <- log(0.2)
    log_sigma_c = as.numeric(glmmres$sdr$par.fixed["theta"])
    
    if (!is.null(sigma_init)) {
      if (!is.na(sigma_init[1])) log_sigma_a <- log(sigma_init[1])
      if (!is.na(sigma_init[2])) log_sigma_b <- log(sigma_init[2])
      if (!is.na(sigma_init[3])) log_sigma_c <- log(sigma_init[3])
    }
    
    trans_rho_a <- trans_rho_b <- trans_rho_c <- 0
    log_phi = log(0.5)
    
    # if (is.null(dat_tuningCPUE)) {
    #   log_q = 0
    #   log_tau = log(0.2)
    # } else {
    #   log_q = mean(log(tuningCPUE))-c_mean
    #   log_tau = log(0.5)
    # }
    
    param = list(a_y=a_y,b_y=b_y,c_y=c_y,a_mean=a_mean,b_mean=b_mean,c_mean=c_mean,
                 log_sigma_a=log_sigma_a,log_sigma_b=log_sigma_b,log_sigma_c=log_sigma_c,
                 trans_rho_a=trans_rho_a,trans_rho_b=trans_rho_b,trans_rho_c=trans_rho_c,
                 log_phi=log_phi,
                 trans_psi=rep(0,3),trans_rho=0
                 #, log_q=log_q,log_tau=log_tau
                 )
  } else {
    param = p0_list
  }
  
  map = list()
  random_c = c()
  if (family=="poisson") map$log_phi = factor(NA)
  
  if (use_MVN==0) {
    if (a_key == 0) { # fixed effect
      map$log_sigma_a = factor(NA)
      map$a_mean = factor(NA)
    } else {
      if (a_key == 1) {
        map$log_sigma_a = factor(NA)
        map$a_y = rep(factor(NA),ny)
      } else {
        if (a_key == 4) map$a_mean = factor(NA)
      }
    }
    if (a_key > 1) random_c = c(random_c,"a_y")
    if (a_key != 3) map$trans_rho_a = factor(NA)
    
    if (b_key == 0) { # fixed effect
      map$log_sigma_b = factor(NA)
      map$b_mean = factor(NA)
    } else {
      if (b_key == 1) {
        map$log_sigma_b = factor(NA)
        map$b_y = rep(factor(NA),ny)
      } else {
        if (b_key == 4) map$b_mean = factor(NA)
      }
    }
    if (b_key > 1) random_c = c(random_c,"b_y")
    if (b_key != 3) map$trans_rho_b = factor(NA)
    
    if (c_key == 0) { # fixed effect
      map$log_sigma_c = factor(NA)
      map$c_mean = factor(NA)
    } else {
      if (c_key == 1) {
        map$log_sigma_c = factor(NA)
        map$c_y = rep(factor(NA),ny)
      } else {
        if (c_key == 4) map$c_mean = factor(NA)
      }
    }
    
    if (c_key > 1) random_c = c(random_c,"c_y")
    if (c_key != 3) map$trans_rho_c = factor(NA)
    
    map$trans_rho = factor(NA)
    map$trans_psi = rep(factor(NA),3)
    
  } else {
    if (use_MVN==1) {
      stop("'use_MVN==1' not allowed (Use 0 or 2-4)")
    } else {
      random_c = c("a_y","b_y","c_y")
      map$trans_rho_a = factor(NA)
      map$trans_rho_b = factor(NA)
      map$trans_rho_c = factor(NA)
      if (use_MVN == 2) { ## white noise
        map$trans_rho = factor(NA)
      } else {
        if (use_MVN == 4) { ## Random walk
          map$a_mean = factor(NA)
          map$b_mean = factor(NA)
          map$c_mean = factor(NA)
          map$trans_rho = factor(NA)
        }
      }
      if (length(map_psi) != 3) stop("The length of 'map_psi' must be three!!!")
      map$trans_psi = factor(map_psi)
    }
  }
  
  if (!is.null(map_add)) {
    tmp = messir::make_named_list(map_add)
    for(i in 1:length(tmp$map_add)) {
      map[[names(tmp$map_add)[i]]] <- map_add[[i]]
    }
  }
  
  obj = TMB::MakeADFun(data=data,parameters=param,random=random_c,
                       DLL="surfnet_peak_model",map=map,silent=silent)
  opt = try(nlminb(obj$par,obj$fn,obj$gr)) 
  if (class(opt) == "try-error") stop("Error in 'nlminb'") 
  if (a_key<2 && b_key<2 && c_key<2 && use_MVN==0) {
    rep = try(sdreport(obj))
  } else {
    rep = try(sdreport(obj,bias.correct = bias_correct))
  }
  if (class(rep) == "try-error") stop("Error in 'sdreport'") 
  
  Res$obj = obj
  Res$opt = opt
  Res$rep = rep
  Res$convergence <- opt$convergence
  Res$pdHess <- rep$pdHess
  Res$data <- data
  Res$param_init <- param
  Res$map <- map
  Res$par_list = obj$env$parList(opt$par)
  
  Res$loglik <- loglik <- -opt$objective
  Res$N <- N <- nrow(dat)
  Res$npar <- npar <- length(opt$par)
  Res$AIC <- AIC <- -2*loglik + 2*npar
  Res$AICc <- AICc <- AIC+2*npar*(npar+1)/(N-npar-1)
  Res$BIC <- BIC <- -2*loglik + npar*log(N)
  
  rep_summary = summary(rep,"report")
  Res$rep_summary = rep_summary
  
    # if (bias_correct) {
  #   VALUE = rep$unbiased$value
  # } else {
    VALUE = rep$value
  # }
  Res$Catch_pred <- Catch_pred <-  VALUE[names(VALUE)=="Catch_pred"]
  Res$CPUE_pred <- CPUE_pred <- VALUE[names(VALUE)=="CPUE_pred"]
  
  dat_pred = dat %>%
    mutate(Catch_pred = Catch_pred, CPUE_pred = CPUE_pred) %>%
    mutate(Catch_se = rep_summary[rownames(rep_summary)=="Catch_pred",2]) %>%
    mutate(CPUE_se = rep_summary[rownames(rep_summary)=="CPUE_pred",2]) %>%
    mutate(Catch_CV = Catch_se/Catch_pred, CPUE_CV = CPUE_se/CPUE_pred) %>%
    mutate(loglik = rep_summary[rownames(rep_summary)=="loglik",1],
           deviance_residual = rep_summary[rownames(rep_summary)=="deviance_residual",1]) %>%
    mutate(Var = rep_summary[rownames(rep_summary)=="Var",1])
  
  Res$dat_pred <- dat_pred
  
  icol = ifelse(bias_correct && any(a_key>1,b_key>1,c_key>1), 3, 1)
  index = data.frame(Year = min(dat$Year):max(dat$Year), 
                     Estimate = rep_summary[rownames(rep_summary)=="index",icol],
                     SE = rep_summary[rownames(rep_summary)=="index",2]) %>%
    dplyr::mutate(CV = SE/Estimate) %>%
    dplyr::mutate(a_par = rep_summary[rownames(rep_summary) == "a_par",1],
                  b_par = rep_summary[rownames(rep_summary) == "b_par",1],
                  c_par = rep_summary[rownames(rep_summary) == "c_par",1])
  Res$index <- index
  
  index_scaled = data.frame(Year = min(dat$Year):max(dat$Year), 
                     Estimate = rep_summary[rownames(rep_summary)=="index_scaled",icol],
                     SE = rep_summary[rownames(rep_summary)=="index_scaled",2]) %>%
    dplyr::mutate(CV = SE/Estimate) %>%
    dplyr::mutate(a_par = rep_summary[rownames(rep_summary) == "a_par",1],
                  b_par = rep_summary[rownames(rep_summary) == "b_par",1],
                  c_par = rep_summary[rownames(rep_summary) == "c_par",1])
  Res$index_scaled <- index_scaled
  
  abc_mean = rep_summary[rownames(rep_summary) == "abc_mean",1]
  names(abc_mean) <- c("a","b","c")
  abc_var = matrix(rep_summary[rownames(rep_summary) == "abc_var",1],ncol=3)
  Res$abc_mean = abc_mean
  Res$abc_var = abc_var
  Res$abc_cor = stats::cov2cor(abc_var)
  
  psi = rep_summary[rownames(rep_summary) == "psi",1]
  names(psi) <- c("a","b","c")
  rho = rep_summary[rownames(rep_summary) == "rho",1]
  Res$psi = psi
  Res$rho = rho

  nTime = max(dat$Time)-min(dat$Time)+1
  CPUE_per_day = data.frame(Year = rep(min(dat$Year):max(dat$Year),nTime),
                            Time = rep(min(dat$Time):max(dat$Time),each=ny)) %>%
    dplyr::mutate(Estimate = rep_summary[rownames(rep_summary)=="cpue_per_day",1],
           SE = rep_summary[rownames(rep_summary)=="cpue_per_day",2])
  if (bias_correct && any(a_key>1,b_key>1,c_key>1)) {
    CPUE_per_day = CPUE_per_day %>% 
      dplyr::mutate(Estimate_unbiased = rep_summary[rownames(rep_summary)=="cpue_per_day",3])
  }
  CPUE_per_day = CPUE_per_day %>% arrange(Year,Time)
  Res$CPUE_per_day <- CPUE_per_day
  
  if (0) {
    g1 = ggplot(data=CPUE_per_day,aes(x=Time)) + 
      geom_area(aes(y=Estimate),fill="red",alpha=0.3)+
      geom_path(aes(y=Estimate,group=Year),size=1,colour="red")+
      geom_point(data = dat_pred, aes(x=Time,y=Catch/Effort,group=Year))+
      facet_wrap(vars(Year),ncol=4, scales="free_y")+
      frasyr::theme_SH()
    g1
    
    # g2 = ggplot(data=CPUE_per_day,aes(x=Time)) + 
    #   geom_area(aes(y=Estimate_unbiased),fill="red",alpha=0.3)+
    #   geom_path(aes(y=Estimate_unbiased,group=Year),size=1,colour="red")+
    #   geom_point(data = dat_pred, aes(x=Time,y=Catch/Effort,group=Year))+
    #   facet_wrap(vars(Year),ncol=4, scales="free_y")+
    #   frasyr::theme_SH()
    # g2
    
    nominal = dat %>% group_by(Year) %>%
      summarise(CPUE = mean(Catch/Effort)) %>%
      mutate(mean = mean(CPUE)) %>%
      mutate(CPUE_scaled = CPUE/mean)
    g3 = ggplot(data=index_scaled,aes(y=Estimate,x=Year)) +
      geom_path() +
      geom_point(data=nominal,aes(x=Year,y=CPUE_scaled))
    g3
  }
  return(Res)
}

# ggsave(g1, file = "CPUE_fitting.png",dpi=600,unit="mm",height=150,width=240)
# ggsave(g1, file = "CPUE_fitting_poisson.png",dpi=600,unit="mm",height=150,width=240)

plot_spm_fit = function(model, CI = 0.95, ncol=4, scales="free_y") {
  CPUE_per_day = model$CPUE_per_day
  dat_pred =model$dat_pred
  CPUE_per_day = CPUE_per_day %>% mutate(CV = SE/Estimate) %>%
    mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>%
    mutate(Lower = Estimate/Cz, Upper = Estimate*Cz)
  g1 = ggplot(data=CPUE_per_day,aes(x=Time)) + 
    geom_ribbon(aes(ymax=Upper,ymin=Lower),fill = "red",alpha=0.3)+
    geom_path(aes(y=Estimate,group=Year),size=1,colour="red")+
    geom_point(data = dat_pred, aes(x=Time,y=Catch/Effort,group=Year))+
    facet_wrap(vars(Year),ncol=ncol, scales=scales)+
    xlab("Days elapsed since May 1")+ylab("CPUE")+
    frasyr::theme_SH()
  return(g1)
}

plot_spm_trend = function(model, CI = 0.95) {
  nominal = model$input$dat %>% group_by(Year) %>%
    summarise(CPUE_mean = mean(Catch/Effort),
              CPUE_max = max(Catch/Effort)) %>%
    mutate(mean_mean = mean(CPUE_mean),
           mean_max = mean(CPUE_max)) %>%
    mutate(Mean = CPUE_mean/mean_mean,
           Max = CPUE_max/mean_max) %>%
    pivot_longer(cols=c(Mean,Max),names_to="Nominal",values_to = "value")
  CI = 0.95
  index_scaled = model$index_scaled %>%
    mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>%
    mutate(Lower = Estimate/Cz, Upper = Estimate*Cz)
  
  g3 = ggplot(data=index_scaled,aes(x=Year)) +
    geom_ribbon(aes(ymax=Upper,ymin=Lower),fill="red",alpha=0.3)+
    geom_path(aes(y=Estimate),size=2,colour="red") +
    geom_point(data=nominal,aes(x=Year,y=value,colour=Nominal,shape=Nominal),size=3)+
    ylab("Scaled Index")+
    frasyr::theme_SH()
  return(g3)
}

plot_spm_trend2 = function(model, CI = 0.95, legend_position = "bottom",
                           legend_direction="horizontal",rec_est=NULL,
                           fill_palette="Set1") {
  nominal = model$input$dat %>% group_by(Year) %>%
    summarise(nominal_mean = mean(Catch/Effort),nominal_max=max(Catch/Effort)) %>%
    dplyr::mutate(Nominal_mean = nominal_mean/mean(nominal_mean),
                  Nominal_max = nominal_max/mean(nominal_max)) %>%
    dplyr::select(Year,Nominal_mean,Nominal_max) %>%
    tidyr::pivot_longer(-Year,names_to="Type",values_to="Scaled_CPUE")
  CI = CI
  index_scaled = model$index_scaled %>%
    mutate(Cz = exp(qnorm(CI+(1-CI)/2)*sqrt(log(1+CV^2)))) %>%
    mutate(Lower = Estimate/Cz, Upper = Estimate*Cz) %>%
    dplyr::select(Year,Estimate,Lower,Upper) %>%
    rename(Scaled_CPUE = Estimate) %>%
    dplyr::mutate(Type = "Standardized")
  
  combined = full_join(nominal,index_scaled) %>%
    dplyr::mutate(Type_f = factor(Type, levels=c("Standardized","Nominal_mean","Nominal_max"))) %>%
    dplyr::select(-Type) %>%
    rename(Type = Type_f)
  
  # labeli = c(Standardized="標準化",Nominal_mean="ノミナル(平均)",Nominal_max="ノミナル(最大)")
  labeli = c(Standardized="Standardized",Nominal_mean="Nominal mean",Nominal_max="Nominal max")
  
  ymax = c(combined$Upper,combined$Scaled_CPUE) %>% max(na.rm=T)*1.1
  g3 = ggplot(data=NULL,aes(x=Year))+ylim(0,ymax)
  
  if (!is.null(rec_est)) {
    trend_mean = combined %>% dplyr::filter(Year %in% rec_est$Year) %>%
      summarise(Average = mean(Scaled_CPUE))
    scale_f = mean(rec_est$Natural)/trend_mean$Average
    rec_est2 = rec_est %>% dplyr::select(-Total) %>% 
      pivot_longer(cols=-Year,names_to="Type",values_to="Recruitment")
    g3 = g3+geom_bar(data=rec_est2,aes(y=Recruitment/scale_f,fill=Type),stat="identity")+
      # scale_fill_brewer(palette="Set3",name="")+
      scale_y_continuous(limits=c(0,ymax),
                         sec.axis = sec_axis(~ . * scale_f,name="Recruitment Estimate"))+
      scale_fill_brewer(palette=fill_palette,name="")
  }
  
  if (CI>0) {
    g3 = g3 + geom_ribbon(data=combined,aes(ymax=Upper,ymin=Lower,fill=Type),alpha=0.4)
  }
  g3 = g3 + geom_path(data=combined,aes(y=Scaled_CPUE,colour=Type,linetype=Type),size=1.5) +
    geom_point(data=combined,aes(y=Scaled_CPUE,colour=Type,shape=Type),size=3) +
    scale_colour_brewer(palette="Set1",name="",labels=labeli)+
    scale_linetype(name="",labels=labeli)+
    scale_shape(name="",labels=labeli)+
    ylab("Scaled Index")+
    # xlab("年")+
    theme_bw(base_size=16)+
    theme(legend.direction=legend_direction,legend.position = legend_position)
  if (is.null(rec_est)) {
    g3 = g3+scale_fill_brewer(palette=fill_palette,guide=FALSE)
  }
  g3
  # +frasyr::theme_SH(legend.position = legend_position)
  invisible(list(data=combined,graph=g3))
  # return(g3)
}


spm_CV = function(model_list, k=10, seed=1,type=0) {
  dat = model_list[[1]]$input$dat
  w_vec = rep_len(1:k, length = nrow(dat))
  set.seed(seed)
  w_vec2 = sample(w_vec)
  jmax = length(model_list)
  cv_summary = NULL
  for (i in 1:k) {
    for (j in 1:jmax) {
      model = model_list[[j]]
      w = ifelse(w_vec2 == i, 0, 1)
      input = model$input
      input$w <- w
      input$p0_list = model$par_list
      res = try(do.call(spm, input))
      if(class(res)=="try-error"){
        rmse = 9999
        mae = 9999
        mad = 9999
        ls = 9999
      } else {
        Var = res$dat_pred$Var[w==0]
        Obs = res$dat_pred$Catch[w==0]
        Pred = res$dat_pred$Catch_pred[w==0]
        CPUE_pred = res$dat_pred$CPUE_pred[w==0]
        CPUE_obs = res$dat_pred$Catch[w==0]/res$dat_pred$Effort[w==0]
        Residual = res$dat_pred$deviance_residual[w==0]
        if (type==0) {
          # standardized by SD
          Error = (Obs-Pred)/sqrt(Var)
        } else {
          if (type==1) {
            #difference of catch 
            Error = (Obs-Pred)
          } else {
            if (type==2) {
              #difference of CPUE
              Error = CPUE_obs-CPUE_pred
            } else {
              # using deviance residuals
              Error = Residual
            }
          }
        }
        nll = -res$dat_pred$loglik[w==0]
        rmse = sqrt(mean(Error^2))
        mae = mean(abs(Error))
        mad = median(abs(Error))
        ls = mean(nll)
      }
      summary = c(model = j, sample_id = i, rmse = rmse, mae = mae, mad=mad, ls = ls)
      cv_summary = rbind(cv_summary, summary)
      print(summary)
    }
  }
  cv_summary_table = cv_summary %>% as_tibble() %>% group_by(model) %>%
    summarise(mean_rmse = mean(rmse), sd_rmse = sd(rmse), 
              mean_mae = mean(mae), sd_mae = sd(mae),
              mean_mad = mean(mad), sd_mad = sd(mad),
              mean_ls = mean(ls), sd_ls = sd(ls))
  return(list(all=cv_summary,summary=cv_summary_table))
}

convert_to_df = function(vpares) {
  Year = as.numeric(colnames(vpares$naa))
  naa = vpares$naa %>% as.data.frame() %>% mutate(Age=0:3) %>% 
    pivot_longer(cols=-Age,names_to="Year",values_to="naa")
  baa = vpares$baa %>% as.data.frame() %>% mutate(Age=0:3) %>% 
    pivot_longer(cols=-Age,names_to="Year",values_to="baa")
  ssb = vpares$ssb %>% as.data.frame() %>% mutate(Age=0:3) %>% 
    pivot_longer(cols=-Age,names_to="Year",values_to="ssb")
  faa = vpares$faa %>% as.data.frame() %>% mutate(Age=0:3) %>% 
    pivot_longer(cols=-Age,names_to="Year",values_to="faa")
  caa = vpares$input$dat$caa %>% as.data.frame() %>% mutate(Age=0:3) %>% 
    pivot_longer(cols=-Age,names_to="Year",values_to="caa")
  df = full_join(naa,baa) %>% full_join(ssb) %>% full_join(faa) %>%
    full_join(caa) %>%
    mutate(Year = as.numeric(Year)) %>%
    arrange(Year,Age)
  
  return(df)
}


# vpares = vout1
# index_value=rev(index_age0$Estimate)[1]
# rec_additive = rec_est0$Additive

one_year_forecast = function(vpares,index_value=NULL,rec_ave_year=5, rec_additive) {
  naa = vpares$naa[,ncol(vpares$naa)]
  faa = vpares$faa[,ncol(vpares$faa)]
  M = vpares$input$dat$M[,ncol(vpares$naa)]
  maa = vpares$input$dat$maa[,ncol(vpares$naa)]
  waa = vpares$input$dat$waa[,ncol(vpares$naa)]
  naa_forecast = c(0,(naa*exp(-M-faa))[-4])+c(0,0,0,(naa*exp(-M-faa))[4])
  if (is.null(index_value)) {
    rec = as.numeric(rowMeans(rev(vpares$naa[1,])[1:rec_ave_year]))
  } else {
    rec_nat = index_value/vpares$q[2]
    # rec_add_dat = data.frame(Rec = unlist(vpares$naa[1,]),Year = as.numeric(colnames(vpares$naa))) %>% full_join(prop_nat) %>%
    #   na.omit() %>% mutate(Rec_additive = Rec*(1-Prop_nat))
    rec_add = mean(rev(rec_additive)[1:rec_ave_year])
    rec = rec_nat + rec_add
  }
  naa_forecast[1] <- rec
  baa_forecast = naa_forecast*waa
  ssb_forecast = baa_forecast*maa
  # res = list(naa=naa_forecast,baa=baa_forecast,ssb=ssb_forecast)
  res = data.frame(Age = 0:3,Year=max(as.numeric(colnames(vpares$naa)))+1,
                   naa=naa_forecast,baa=baa_forecast,ssb=ssb_forecast,faa=vpares$Fc.at.age) %>%
    dplyr::mutate(caa=naa*exp(-0.5*M)*(1-exp(-faa))) %>%
    as_tibble()
  vpa_df = convert_to_df(vpares)
  Res = full_join(vpa_df,res)
  
  return(Res)
}

# dat_forecast = vout1_forecast
# dat_forecast = vout_forecast
get_mohn_graph = function(dat_forecast,
                          what_plot=c("Recruit","Biomass","SSB","F","U","Catch"),
                          ncol=3, number_unit=1000,biomass_unit=1000,
                          plot_start_year = 2005) {
  dat_forecast = dat_forecast %>% 
    dplyr::mutate(retro_year = factor(retro_year),cwaa = caa*baa/naa)
  dat_r = dat_forecast %>% dplyr::filter(Age==0) %>%
    dplyr::select(Year,naa,retro_year) %>%
    rename(Recruit = naa) %>% mutate(Recruit=Recruit/number_unit)
  dat_n = dat_forecast %>% group_by(Year,retro_year) %>%
    summarise(Number = sum(naa)) %>% mutate(Number=Number/number_unit)
  dat_b = dat_forecast %>% group_by(Year,retro_year) %>%
    summarise(Biomass = sum(baa)) %>% mutate(Biomass=Biomass/biomass_unit)
  dat_s = dat_forecast %>% group_by(Year,retro_year) %>%
    summarise(SSB = sum(ssb)) %>% mutate(SSB=SSB/biomass_unit)
  dat_f = dat_forecast %>% group_by(Year,retro_year) %>%
    summarise(F = mean(faa))
  dat_u = dat_forecast %>% group_by(Year,retro_year) %>%
    summarise(U = sum(cwaa)/sum(baa))
  dat_c = dat_forecast %>% group_by(Year,retro_year) %>%
    summarise(Catch = sum(cwaa)) %>% mutate(Catch=Catch/biomass_unit)
  dat_combined = full_join(dat_n,dat_b) %>% full_join(dat_r) %>% full_join(dat_s) %>% 
    full_join(dat_f) %>% full_join(dat_u) %>% full_join(dat_c)
  # View(dat_combined)
  dat_longer = dat_combined %>% 
    pivot_longer(cols=c(-Year,-retro_year),names_to="Variable",values_to="Value") %>%
    dplyr::mutate(Variable=factor(Variable,levels=c("Recruit","Number","Biomass","SSB","F","U","Catch"))) %>%
    arrange(retro_year,Variable,Year) %>% dplyr::filter(Variable %in% what_plot)
  dat_plot = filter(dat_longer,Year>plot_start_year)
  dat_ymax = dat_plot %>% group_by(Variable) %>% summarise(ymax=max(Value)*1.15)
  g1 = ggplot(data=dat_plot)+
    geom_path(aes(x=Year,y=Value,group=retro_year,colour=retro_year,linetype=retro_year),show.legend=FALSE,size=1)+
    geom_point(aes(x=Year,y=Value,group=retro_year,colour=retro_year),size=1.5,show.legend=FALSE)+
    geom_blank(data=dat_ymax,
               aes(y=ymax))+
    facet_wrap(vars(Variable),scales="free_y",ncol=ncol)+
    theme_bw(base_size=16)+ylab("")
  # g1
  
  max_year = dat_combined$Year %>% max()
  n_retro = dat_combined$retro_year %>% as.character() %>% as.numeric() %>% max()
  rho_list = sapply(1:n_retro, function(i) {
    base = dat_combined %>% filter(retro_year == 0 & Year == max_year-i) %>%
      ungroup() %>% 
      dplyr::select(-Year,-retro_year) %>% unlist()
    est = dat_combined %>% filter(retro_year == i & Year == max_year-i) %>%
      ungroup() %>% 
      dplyr::select(-Year,-retro_year) %>% unlist()
    rho = est/base-1
    rho
  })
  rho=rowMeans(as.matrix(rho_list))
  rho_data = data.frame(Variable=names(rho),rho=rho) %>% 
    dplyr::mutate(Variable=factor(Variable,levels=c("Recruit","Number","Biomass","SSB","F","U","Catch"))) %>%
    dplyr::filter(Variable %in% what_plot) %>%
    full_join(dat_ymax) %>% 
    rename(Value=ymax) %>%
    mutate(Year = max_year) %>%
    mutate(label=sprintf("rho == %.2f",rho))
  g1 = g1 + geom_text(data=rho_data,parse=TRUE,aes(x=Year,y=Value,label=label,hjust=1,vjust=1))+
    ylim(0,NA)
  g1
  
  RES = list(dat_combined=dat_combined,dat_longer=dat_longer,dat_plot=dat_plot,gg=g1,rho_list=rho_list,rho=rho)
}


