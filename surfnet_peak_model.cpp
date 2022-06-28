
#include <TMB.hpp>
#include <iostream>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(Catch);  // size: sample size of data
  DATA_VECTOR(Effort);  // size: sample size of data
  DATA_IVECTOR(iYear);  // size: sample size of data
  DATA_VECTOR(Time);  // size: sample size of data
  DATA_INTEGER(a_key);   //0: fixed effect, 1: constant, 2: white noise, 3: AR1, 4: Random walk
  DATA_INTEGER(b_key);
  DATA_INTEGER(c_key);
  DATA_INTEGER(family); //0: Poisson, 1: negative binomial(variance = mu(1+phi)), 2: negative binomial (variance = mu(1+phi*mu))
  DATA_VECTOR(Time_grid);
  DATA_IVECTOR(w);  // weight of each data (usually all 1)
  DATA_SCALAR(Time_sd);
  DATA_INTEGER(use_MVN); //whether using MVN or not,  0: not use, 1 not recognized, 2: white noise, 3: multivariate AR1, 4: Random walk

  // DATA_INTEGER(tuningCPUE_key); // whether to use Kosoko CPUE for tuning
  // DATA_VECTOR(tuningCPUE);
  // DATA_IVECTOR(iYear_tuningCPUE);

  PARAMETER_VECTOR(a_y);   // size: length of years (random effect)
  PARAMETER_VECTOR(b_y);   // size: length of years (random effect)
  PARAMETER_VECTOR(c_y);   // size: length of years (random effect)
  PARAMETER(a_mean);
  PARAMETER(b_mean);
  PARAMETER(c_mean);
  PARAMETER(log_sigma_a);
  PARAMETER(log_sigma_b);
  PARAMETER(log_sigma_c);
  PARAMETER(trans_rho_a);
  PARAMETER(trans_rho_b);
  PARAMETER(trans_rho_c);
  PARAMETER(log_phi);   // overdispersion factor for a negative binomial distribution
  PARAMETER_VECTOR(trans_psi); //correlation coefficients among a_y, b_y, and c_y (length is 3: first=a_y&b_y, second=b_y&c_y, third=c_y&a_y)
  PARAMETER(trans_rho); //correlation coefficients when use_MVN=3 (multivariate AR1)

  // PARAMETER(log_q);
  // PARAMETER(log_tau);

  Type sigma_a = exp(log_sigma_a);
  Type sigma_b = exp(log_sigma_b);
  Type sigma_c = exp(log_sigma_c);

  Type rho_a = (exp(trans_rho_a)-Type(1.0))/(exp(trans_rho_a)+Type(1.0));
  Type rho_b = (exp(trans_rho_b)-Type(1.0))/(exp(trans_rho_b)+Type(1.0));
  Type rho_c = (exp(trans_rho_c)-Type(1.0))/(exp(trans_rho_c)+Type(1.0));

  Type rho = (exp(trans_rho)-Type(1.0))/(exp(trans_rho)+Type(1.0));
  vector<Type> psi = (exp(trans_psi)-Type(1.0))/(exp(trans_psi)+Type(1.0));

  Type phi = exp(log_phi);

  vector<Type> a_par(c_y.size());
  vector<Type> b_par(c_y.size());
  vector<Type> c_par(c_y.size());

  array<Type> abc_mat(3,c_y.size());
  for (int i=0;i<c_y.size();i++) {
    abc_mat(0,i) = a_y(i);
    abc_mat(1,i) = b_y(i);
    abc_mat(2,i) = c_y(i);
  }

  vector<Type> abc_mean(3);
  abc_mean(0) = a_mean;
  abc_mean(1) = b_mean;
  abc_mean(2) = c_mean;

  array<Type> abc_diff(3,c_y.size());
  for (int i=0;i<c_y.size();i++) {
    abc_diff.col(i) = abc_mat(i)-abc_mean;
  }

  matrix<Type> abc_var(3,3);
  abc_var(0,0) = pow(sigma_a,Type(2.0)); //variance of a_y
  abc_var(1,1) = pow(sigma_b,Type(2.0)); //variance of b_y
  abc_var(2,2) = pow(sigma_c,Type(2.0)); //variance of b_y
  abc_var(0,1) = sigma_a*sigma_b*psi(0); //covariance of a_y and b_y
  abc_var(1,0) = sigma_a*sigma_b*psi(0); //covariance of a_y and b_y
  abc_var(1,2) = sigma_b*sigma_c*psi(1); //covariance of a_y and b_y
  abc_var(2,1) = sigma_b*sigma_c*psi(1); //covariance of a_y and b_y
  abc_var(0,2) = sigma_a*sigma_c*psi(2); //covariance of a_y and b_y
  abc_var(2,0) = sigma_a*sigma_c*psi(2); //covariance of a_y and b_y

  using namespace density;

  MVNORM_t<Type> neg_log_density_abc(abc_var);
  Type nll = 0.0;

  if (use_MVN == 0) {
    // parameterization of a
    // process likelihood
    if (a_key == 2) { // white noise
      nll -= sum(dnorm(a_y, a_mean, sigma_a, true));
    }
    if (a_key == 3) { // AR1
      nll += SCALE(AR1(rho_a),sigma_a)(a_y-a_mean);
    }
    if (a_key == 4) { // random walk
      for (int i=1;i<c_y.size();i++) {
        nll -= dnorm(a_y(i), a_y(i-1), sigma_a, true);
      }
    }

    if (a_key == 1) { //constant
      for (int i=0;i<c_y.size();i++) {
        a_par(i) = a_mean;
      }
    } else {
      a_par = a_y;
    }

    // parameterization of b
    // process likelihood
    if (b_key == 2) { // white noise
      nll -= sum(dnorm(b_y, b_mean, sigma_b, true));
    }
    if (b_key == 3) { // AR1
      nll += SCALE(AR1(rho_b),sigma_b)(b_y-b_mean);
    }
    if (b_key == 4) { // random walk
      for (int i=1;i<c_y.size();i++) {
        nll -= dnorm(b_y(i), b_y(i-1), sigma_b, true);
      }
    }

    if (b_key == 1) { //constant
      for (int i=0;i<c_y.size();i++) {
        b_par(i) = b_mean;
      }
    } else {
      b_par = b_y;
    }

    // parameterization of c
    // process likelihood
    if (c_key == 2) { // white noise
      nll -= sum(dnorm(c_y, c_mean, sigma_c, true));
    }
    if (c_key == 3) { // AR1
      nll += SCALE(AR1(rho_c),sigma_c)(c_y-c_mean);
    }
    if (c_key == 4) { // random walk
      for (int i=1;i<c_y.size();i++) {
        nll -= dnorm(c_y(i), c_y(i-1), sigma_c, true);
      }
    }

    if (c_key == 1) { //constant
      for (int i=0;i<c_y.size();i++) {
        c_par(i) = c_mean;
      }
    } else {
      c_par = c_y;
    }
  } else {
    if (use_MVN == 1) {
      error("use_MVN code not recognized");
    } else {
      if (use_MVN == 2) { //white noise
        for (int i=0;i<c_y.size();i++) {
          nll += neg_log_density_abc(abc_diff.col(i));
        }
      } else{
        if (use_MVN == 3) { //AR1 process
          nll += AR1(rho,neg_log_density_abc)(abc_diff);
        } else {
          for (int i=1;i<c_y.size();i++) {
            nll += neg_log_density_abc(abc_mat.col(i)-abc_mat.col(i-1));
          }
        }
      }
    }
    a_par = a_y;
    b_par = b_y;
    c_par = c_y;
  }


  // observation likelihood
  vector<Type> log_CPUE_pred(Catch.size());
  vector<Type> CPUE_pred(Catch.size());
  vector<Type> log_Catch_pred(Catch.size());
  vector<Type> Catch_pred(Catch.size());
  vector<Type> Var(Catch.size());
  vector<Type> Var_obs(Catch.size());
  vector<Type> loglik(Catch.size());
  vector<Type> deviance_residual(Catch.size());

  for (int i=0;i<Catch.size();i++) {
    log_CPUE_pred(i) = -exp(a_par(iYear(i)))*pow(Time(i)-b_par(iYear(i)),Type(2.0))+c_par(iYear(i));
    CPUE_pred(i) = exp(log_CPUE_pred(i));
    log_Catch_pred(i) = log_CPUE_pred(i) + log(Effort(i));
    Catch_pred(i) = exp(log_Catch_pred(i));
    if (family == 0) { // Poisson distribution
      Var(i) = Catch_pred(i);
      Var_obs(i) = Catch(i);
      loglik(i) = dpois(Catch(i),Catch_pred(i),true);
    if (Catch(i)>0.0) {
      deviance_residual(i) = Type(2.0)*(dpois(Catch(i),Catch(i),true)-loglik(i));
      SIMULATE{
        Catch(i) = rpois(Catch_pred(i));
      }
    } else {
      deviance_residual(i) = Type(2.0)*(-loglik(i));
    }
    } else { if (family == 1) { //nbinom1
      Var(i) = Catch_pred(i)*(Type(1.0)+phi);
      Var_obs(i) = Catch(i)*(Type(1.0)+phi);
      loglik(i) = dnbinom_robust(Catch(i),log_Catch_pred(i),log_phi+log_Catch_pred(i),true);
    } else {
      Var(i) = Catch_pred(i)*(Type(1.0)+phi*Catch_pred(i));
      Var_obs(i) = Catch(i)*(Type(1.0)+phi*Catch(i));
      loglik(i) = dnbinom_robust(Catch(i),log_Catch_pred(i),log_phi+Type(2.0)*log_Catch_pred(i),true);
    }
    // loglik(i) = dnbinom2(Catch(i),Catch_pred(i),Var(i),true);
    // loglik(i) = dnbinom(Catch(i),pow(Catch_pred(i),Type(2.0))/(Var(i)-Catch_pred(i)),Catch_pred(i)/Var(i),true);
    if (Catch(i)>0.0) {
      // deviance_residual(i) = Type(2.0)*(dnbinom2(Catch(i),Catch(i),Var(i),true)-loglik(i));
      deviance_residual(i) = dnbinom(Catch(i),pow(Catch_pred(i),Type(2.0))/(Var(i)-Catch_pred(i)),pow((Catch(i)*(Var(i)-Catch_pred(i))/pow(Catch_pred(i),Type(2.0)))+Type(1.0),Type(-1.0)),true);
      deviance_residual(i) -= loglik(i);
      deviance_residual(i) *= Type(2.0);      //
      // deviance_residual(i) = Type(2.0)*(dnbinom(Catch(i),pow(Catch(i),Type(2.0))/(Var(i)-Catch(i)),Catch(i)/Var(i),true)-loglik(i));
    } else {
      deviance_residual(i) = Type(2.0)*(-loglik(i));
    }
    SIMULATE{
      Catch(i) = rnbinom2(Catch_pred(i),Var(i));
    }
    }
    deviance_residual(i) = pow(deviance_residual(i), Type(0.5));
    deviance_residual(i) *= CppAD::CondExpLe(Catch(i),Catch_pred(i),Type(-1.0),Type(1.0));
    nll -= w(i)*loglik(i);
  }

  // calculation of abundance index
  vector<Type> index(c_y.size());
  //index.fill(0.0);
  array<Type>  cpue_per_day(c_y.size(),Time_grid.size());
  // Type PI=3.14159265358979323846;

  for (int i=0;i<c_y.size();i++) {
    for (int j=0;j<Time_grid.size();j++) {
      cpue_per_day(i,j) = -exp(a_par(i))*pow(Time_grid(j)-b_par(i),Type(2.0))+c_par(i);
      cpue_per_day(i,j) = exp(cpue_per_day(i,j));
  //    index(i) += cpue_per_day(i,j);
    }
     //index(i) /= Time_grid.size(); // calculate average CPUE rather than sum
     index(i) = exp(c_par(i))*pow(M_PI/exp(a_par(i)),Type(0.5));
     index(i) *= Time_sd;
  }

  vector<Type> index_scaled = index/(sum(index)/index.size());

  // vector<Type> loglik_tuning(tuningCPUE.size());
  // loglik_tuning.fill(0.0);
  // vector<Type> deviance_residual_tuning(tuningCPUE.size());
  // deviance_residual_tuning.fill(0.0);
  // Type tau = exp(log_tau);
  //
  // if (tuningCPUE_key==1) {
  //   for (int i=0;i<tuningCPUE.size();i++) {
  //     loglik_tuning(i) += dnorm(log(tuningCPUE(i)),log_q+c_par(iYear_tuningCPUE(i)),tau,true);
  //     deviance_residual_tuning(i) += Type(2.0)*(dnorm(log(tuningCPUE(i)),log(tuningCPUE(i)),tau,true)-loglik_tuning(i));
  //     deviance_residual_tuning(i) = pow(deviance_residual_tuning(i), Type(0.5));
  //     deviance_residual_tuning(i) *= CppAD::CondExpLe(tuningCPUE(i),log_q+c_par(iYear_tuningCPUE(i)),Type(-1.0),Type(1.0));
  //     nll -= loglik_tuning(i);
  //     SIMULATE{
  //       tuningCPUE(i) = exp
  //     }
  //   }
  // }

  // Type process_likelihood = -nll-sum(loglik)-sum(loglik_tuning);
  Type process_likelihood = -nll-sum(loglik);

  SIMULATE {
    REPORT(Catch);
    // REPORT(tuningCPUE);
  }
  ADREPORT(sigma_a);
  ADREPORT(sigma_b);
  ADREPORT(sigma_c);
  ADREPORT(rho_a);
  ADREPORT(rho_b);
  ADREPORT(rho_c);
  ADREPORT(phi);
  ADREPORT(a_par);
  ADREPORT(b_par);
  ADREPORT(c_par);
  ADREPORT(exp(a_par));
  ADREPORT(log_CPUE_pred);
  ADREPORT(CPUE_pred);
  ADREPORT(log_Catch_pred);
  ADREPORT(Catch_pred);
  ADREPORT(Var);
  ADREPORT(Var_obs);
  ADREPORT(loglik);
  ADREPORT(deviance_residual);
  // ADREPORT(tau);
  // ADREPORT(loglik_tuning);
  // ADREPORT(deviance_residual_tuning);
  ADREPORT(process_likelihood);
  ADREPORT(index);
  ADREPORT(index_scaled);
  ADREPORT(cpue_per_day);
  ADREPORT(psi);
  ADREPORT(rho);
  ADREPORT(abc_mat);
  ADREPORT(abc_mean);
  ADREPORT(abc_var);

  return nll;
}
