#---------------------------------------
# Divergent thermal responses of energy acquisition and mobilization reproduce 
#  growth–metabolism decoupling in ectotherms
# Irene Moro-Martinez, Andrea Campos-Candela, Joaquim Tomas-Ferrer & Miquel Palmer
# Submitted to Ecological Modelling
#---------------------------------------
# Last update: 6 nov 2026

#-----
# 1: Loading data and libraries
#-----
remove(list=ls())
library(cmdstanr)
load("input.RData")

#-----       
# 2: STAN model
#-----
sink("DEB.stan")
cat(" // first line
functions {

  ////// DEB Function ////// 
  vector DEB(
    real t, vector y,
    // core parameters
    real f, real E_G, real kap, real p_Am, real p_M, real k_J, real v,
    // temperature parameters
    real Kelvin, real Temp_mean, real Temp_amp, real pi2f, real Temp_phi,
    real T_A, real T_A2, real T_AL, real T_L, real T_AH, real T_H, real dH, real T_ref,
    real birth_correction
  ){
    // --- Derivatives ---  
    vector[2] dydt;  // 1: dEdt; 2: dVdt;

    // --- Precomputing repeated calculations ---
    real V23 = pow(y[2], 2.0/3.0); // surface

    // --- Temperature correction and fluxes ---
    real Temp = Kelvin + Temp_mean + Temp_amp * sin(pi2f * (t + birth_correction) + Temp_phi);
    
      // Model #1
      /*
      // all fluxes: exponential temperature correction
      real cT2 = exp(T_A2/T_ref - T_A2/Temp);
      // --- Fluxes ---
      real pA = f * (p_Am*cT2) * V23;                                                    // assimilation
      real pS = (p_M*cT2) * y[2];                                                        // maintenance
      real pC = (y[1]/y[2]) * (((v*cT2) * E_G * V23 + pS) / (E_G + kap * y[1]/y[2]));    // mobilization
      real pG = (kap*pC - pS);                                                           // growth
      */
    
      // Model #2
      /*
      // all fluxes: unimodal temperature correction
      real sT1 = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/T_H - T_AH/Temp);
      real sT_ref = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref);
      real cT1 = (sT_ref / sT1) * exp(T_A/T_ref - T_A/Temp);
      // --- Fluxes ---
      real pA = f * (p_Am*cT1) * V23;                                                    // assimilation
      real pS = (p_M*cT1) * y[2];                                                        // maintenance
      real pC = (y[1]/y[2]) * (((v*cT1) * E_G * V23 + pS) / (E_G + kap * y[1]/y[2]));    // mobilization
      real pG = (kap*pC - pS);                                                           // growth
      */
    
      // Model #3
      // Assimiation: unimodal temperature correction
      real sT1 = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/T_H - T_AH/Temp);
      real sT_ref = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref);
      real cT1 = (sT_ref / sT1) * exp(T_A/T_ref - T_A/Temp);
      // other fluxes: exponential temperature correction
      real cT2 = exp(T_A2/T_ref - T_A2/Temp);
      // --- Fluxes ---
      real pA = f * (p_Am*cT1) * V23;                                                    // assimilation
      real pS = (p_M*cT2) * y[2];                                                        // maintenance
      real pC = (y[1]/y[2]) * (((v*cT2) * E_G * V23 + pS) / (E_G + kap * y[1]/y[2]));    // mobilization
      real pG = (kap*pC - pS);                                                           // growth
      
      // Model #4
      /*
      // Assimiation: unimodal temperature correction
      real sT1a = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/T_H - T_AH/Temp);
      real sT_refa = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref);
      real cT1 = (sT_refa / sT1a) * exp(T_A/T_ref - T_A/Temp);
      // other fluxes: unimodal temperature correction with differnent T_H
      real sT1b = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/(T_H+dH) - T_AH/Temp);
      real sT_refb = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/(T_H+dH) - T_AH/T_ref);
      real cT2 = (sT_refb / sT1b) * exp(T_A/T_ref - T_A/Temp);;
      // --- Fluxes ---
      real pA = f * (p_Am*cT1) * V23;                                                    // assimilation
      real pS = (p_M*cT2) * y[2];                                                        // maintenance
      real pC = (y[1]/y[2]) * (((v*cT2) * E_G * V23 + pS) / (E_G + kap * y[1]/y[2]));    // mobilization
      real pG = (kap*pC - pS);                                                           // growth 
      */ 
    
    // derivatives
    dydt[1] = pA-pC;
    dydt[2] = pG/E_G;
    
    return dydt;
  } // end DEB function
  
  ////// Function O2_rate ////// 
  array [] real O2_rate(
    // stete variables
    array [] vector y, //E,V and E_H
    // number of samples
    int n,
    // core parameters
    real f, real E_G, real kap, real p_Am, real p_M,real k_J, real v, real kap_R,
    // temperature parameters
    real Kelvin, real Temp_mean, real Temp_amp, real pi2f, real Temp_phi,
    real T_A, real T_A2, real T_AL, real T_L, real T_AH, real T_H, real dH, real T_ref,
    real birth_correction,
    real Temperature_O2,
    // age
    array [] real age,
    // O2 auxiliary parameters
    array [,] real eta_O,array [,] real n_M_n_O

    ){
    
      // expected rate (O2/day)
      array [n] real O2_hat; 

      // Auxiliary variables
      array [3,n] real p_ADG;
      array [4,n] real J_M;
      
      for (i in 1:n){
      
        // --- Precomputing repeated calculations ---
        real V23 = pow(y[i,2], 2.0/3.0); // surface
        
        // --- Temperature correction ---
        real Temp = Kelvin + Temperature_O2;
        
          // Model #1
          /*
          // all fluxes: exponential temperature correction
          real cT2 = exp(T_A2/T_ref - T_A2/Temp);
          // --- Fluxes ---
          real pA = f * (p_Am*cT2) * V23;                                                            // assimilation
          real pS = (p_M*cT2) * y[i,2];                                                              // maintenance
          real pC = (y[i,1]/y[i,2]) * (((v*cT2) * E_G * V23 + pS) / (E_G + kap * y[i,1]/y[i,2]));    // mobilization
          real pG = (kap*pC - pS);                                                                   // growth
          real pD = pS + (1-kap)*pC;                                                                 // dissipation
          */
        
          // Model #2
          /*
          // all fluxes: unimodal temperature correction
          real sT1 = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/T_H - T_AH/Temp);
          real sT_ref = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref);
          real cT1 = (sT_ref / sT1) * exp(T_A/T_ref - T_A/Temp);
          // --- Fluxes ---
          real pA = f * (p_Am*cT1) * V23;                                                            // assimilation
          real pS = (p_M*cT1) * y[i,2];                                                              // maintenance
          real pC = (y[i,1]/y[i,2]) * (((v*cT1) * E_G * V23 + pS) / (E_G + kap * y[i,1]/y[i,2]));    // mobilization
          real pG = (kap*pC - pS);                                                                   // growth
          real pD = pS + (1-kap)*pC;                                                                 // dissipation
          */
        
          // Model #3
          // Assimiation: unimodal temperature correction
          real sT1 = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/T_H - T_AH/Temp);
          real sT_ref = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref);
          real cT1 = (sT_ref / sT1) * exp(T_A/T_ref - T_A/Temp);
          // other fluxes: exponential temperature correction
          real cT2 = exp(T_A2/T_ref - T_A2/Temp);
          // --- Fluxes ---
          real pA = f * (p_Am*cT1) * V23;                                                            // assimilation
          real pS = (p_M*cT2) * y[i,2];                                                              // maintenance
          real pC = (y[i,1]/y[i,2]) * (((v*cT2) * E_G * V23 + pS) / (E_G + kap * y[i,1]/y[i,2]));    // mobilization
          real pG = (kap*pC - pS);                                                                   // growth
          real pD = pS + (1-kap)*pC;                                                                 // dissipation
          
          // Model #4
          /*
          // Assimiation: unimodal temperature correction
          real sT1a = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/T_H - T_AH/Temp);
          real sT_refa = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref);
          real cT1 = (sT_refa / sT1a) * exp(T_A/T_ref - T_A/Temp);
          // other fluxes: unimodal temperature correction with differnent T_H
          real sT1b = 1 + exp(T_AL/Temp - T_AL/T_L) + exp(T_AH/(T_H+dH) - T_AH/Temp);
          real sT_refb = 1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/(T_H+dH) - T_AH/T_ref);
          real cT2 = (sT_refb / sT1b) * exp(T_A/T_ref - T_A/Temp);;
          // --- Fluxes ---
          real pA = f * (p_Am*cT1) * V23;                                                            // assimilation
          real pS = (p_M*cT2) * y[i,2];                                                              // maintenance
          real pC = (y[i,1]/y[i,2]) * (((v*cT2) * E_G * V23 + pS) / (E_G + kap * y[i,1]/y[i,2]));    // mobilization
          real pG = (kap*pC - pS);                                                                   // growth
          real pD = pS + (1-kap)*pC;                                                                 // dissipation
          */ 
          
        // Oxygen 
        p_ADG[1,i] = pA;
        p_ADG[2,i] = pD;
        p_ADG[3,i] = pG;
      }
      
      J_M = to_array_2d(to_matrix(n_M_n_O[1:4,1:4]) * (to_matrix(eta_O[1:4,1:3]) * to_matrix(p_ADG[1:3,1:n])));
      
      for (i in 1:n){   
        O2_hat[i] = 32.0 * J_M[3,i]; // (gr/day) O2 consumption per day
      }
    
      return O2_hat;
    } // end O2 function
  
    ////// Function partial_sum_lpdf (for multicore paralelization) ////// 
    real partial_sum_lpdf(array [,,] real slice_obs,  // 0bservations
                        int start, int end,         // reduce_sum  paramters,
                        
                        //array [] int is_train_fish, // 1/0 for filtering training fish  
                        
                        // arguments for ode_rk45 (DEB model)
                        // state variables at t0, observation dates, number of obsrvations per fish, number of fish
                        array [] vector y0, array [] real t0, array [,] real age, int n, int N, 
                        
                        // DEB main parameters
                        real f, real E_G, real kap, array [] real p_Am, real p_M, real k_J, real v, real kap_R,
                        
                        // temperature-related paramters
                        real Kelvin, array [] real Temp_mean , array [] real Temp_amp , real pi2f, array [] real Temp_phi,
                        real T_A, real T_A2, real T_AL, real T_L, real T_AH, real T_H, real dH, real T_ref,
                        array [] real birth_correction, array [] int ID_temp,
                        array [] real Temperature_O2,
                            
                        // linking state variables with obervable variables
                        real del_M, real d_V, real w_E, real mu_E,
                        array [,] real eta_O, array [,] real n_M_n_O, 
                        
                        // error
                        real sd_length, real sd_weight, real sd_O2

                        ) {

    array [n] vector[2] y;    // numerical integration results
    real likelihood=0;        // cummulated likelihood across variables and fish
    array [n] real O2_hat;

    for (i in 1:end-start+1){  // (fish i; start-to-end defines the especific fish managed by a given core)
        
         y = ode_rk45(DEB, y0[start+i-1], t0[start+i-1] , age[start+i-1,1:n],
                    // DEB main parameters
                    f,E_G,kap,p_Am[start+i-1],p_M,k_J,v,
                    // temperature-related paramters
                    Kelvin,Temp_mean[ID_temp[start+i-1]],Temp_amp[ID_temp[start+i-1]],pi2f,Temp_phi[ID_temp[start+i-1]],
                    T_A,T_A2,T_AL,T_L,T_AH,T_H,dH,T_ref,birth_correction[start+i-1]
        );
       
        O2_hat=O2_rate(y[1:n], // state variables
                      n,
                      // core parameters
                      f,E_G,kap,p_Am[start+i-1],p_M,k_J,v, kap_R,
                      // temperature parameters
                      Kelvin,Temp_mean[ID_temp[start+i-1]],Temp_amp[ID_temp[start+i-1]],pi2f,Temp_phi[ID_temp[start+i-1]],
                      T_A,T_A2,T_AL,T_L,T_AH,T_H,dH,T_ref,birth_correction[start+i-1],
                      Temperature_O2[start+i-1],
                      // age
                      age[start+i-1,1:n],
                      // auxiliary parameters
                      eta_O[,],n_M_n_O[,]
        );
       
        for (j in 1:n){
          if(slice_obs[i,1,j]>0.0){  // missing values as -999
            likelihood += normal_lupdf(slice_obs[i,1,j] | y[j,2]^(1.0/3.0)/del_M , sd_length);
          }
          if(slice_obs[i,2,j]>0.0){
            //likelihood += normal_lupdf(slice_obs[i,2,j] | y[j,2] + w_E/mu_E/d_V*y[j,1] , sd_weight);
            likelihood += lognormal_lupdf(slice_obs[i,2,j] | log(y[j,2] + w_E/mu_E/d_V*y[j,1]) , sd_weight);
          }
          if(slice_obs[i,3,j]>0.0){
            //likelihood += normal_lupdf(slice_obs[i,3,j] | O2_hat[j] , sd_O2);
            likelihood += gamma_lupdf(slice_obs[i,3,j] | square(O2_hat[j])/square(sd_O2), O2_hat[j]/square(sd_O2) );
          }
       }

    } // end fish
    return likelihood;
  }// end reduce sum function
}

data {
  int n;                                               // number of observations per fish
  int N;                                               // number of fish
  array [N,3,n] real obs;                              // observations; fishes,c(length [cm],weight [gr],O2 [gr/day] ), repeated measures (1:n)
  array [N,3] real obs0;                               // observations at t0; fishes,c(length,weight,O2)
  array [N] real t0;                                   // initial time [day]
  array [N,n] real age;                                // times at which the system is observed [day]
  int n_temp;                                          // number of temperature groups
  array [N] int ID_temp;                               // temperaure group membership
  array [n_temp] int n_group;                          // number of fish in each group
  
  // DEB parameters asumed be be known (not estimated) 
  real f;                            // [-] scaled functional response  
  real E_G;                          // [J/cm^3] spec cost for structure
  real kap;                          // [-] allocation fraction to soma
  real p_M;                          // [J/cm3/d] somatic maintenance
  real k_J;                          // [day-1] maturity maint rate coefficient
  real kap_R;                        // [-] reproduction efficiency
  
  // other parameters
  real del_M;                        // [-] shape coefficicent      
  real w_E;                          // [gr/mol] dry weight     
  real mu_E;                         // [J/mol] chemical potential
  real d_V;                          // water content
  array [4,4] real n_M_n_O;          // energy-mass couplers
  array [4,3] real eta_O;            // energy-mass couplers
  
  // temperature parameters
  real Kelvin;                        // Celsius to Kelvin
  array [n_temp] real Temp_mean;      // Sinusoidal function parameter for mean temperature
  array [n_temp] real Temp_amp;       // Sinusoidal function parameter for amplitude
  real pi2f;                          // 2*pi
  array [n_temp] real Temp_phi;       // Sinusoidal function parameter for phase
  real T_ref;                         // [K] reference temperature
  array [N] real birth_correction;    // days from jan 1th
  array [N] real Temperature_O2;      // [C] Temperature at respirometries 
  
  // priors
  array [N] real prior_E0_i;          // Energy at t0
  array [N] real prior_V0_i;          // Volume at t0
  array [2] real prior_p_Am;          // mean for maximum assimulation rate 
  array [2] real prior_p_Am_sd;       // sd for maximum assimulation rate 
  array [2] real prior_v;             // conductance
  array [2] real prior_sd_length;     // length error
  array [2] real prior_sd_weight;     // weight error
  array [2] real prior_sd_O2;         // O2 error
  array [2] real prior_T_A;           // Arrhenius temperature
  array [2] real prior_T_H;           // upper boundary tolerance range
  array [2] real prior_T_AH;          // Arrhenius temperature for upper boundary
  array [2] real prior_T_L;           // lower boundary tolerance range
  array [2] real prior_T_AL;          // Arrhenius temperature for lower boundary
  
  int grainsize;                      // parameter for optimizing multicore parametrization
}

transformed data {
}

parameters {  
  // DEB parameters
  real <lower=0> p_Am_mu;             // [J/cm2/day] mean for maximum assimulation rate 
  real <lower=0> p_Am_sd;             // [J/cm2/day] sd for maximum assimulation rate 
  array [N] real p_Am_star;           // standardized pAm at fish level 
  real <lower=0> v;                   // [cm/day] energy conductance

  // initial values for state variables
  array [N]  real <lower=0> E0_i;     // [J] Energy at t0
  array [N]  real <lower=0> V0_i;     // [cm3] Volume at t0

  // temperature
  real <lower=0.0> T_A;                // [K] Arrhenius temperature for assimilation         
  real <lower=0.0> T_A2;               // [K] Arrhenius temperature for the other fluxes
  real <lower=0.0> T_H;                // [K] upper boundary tolerance range
  real <lower=0.0> T_AH;               // [K] Arrhenius temperature for upper boundary
  real <lower=0.0> T_L;                // [K] lower boundary tolerance range
  real <lower=0.0> T_AL;               // [K] Arrhenius temperature for lower boundary
  real dH;                             // [K] differnece in T_H (Model #4)                          
  
  // observation error
  real <lower=0> sd_length;            // [cm] length
  real <lower=0> sd_weight;            // [gr] weight
  real <lower=0> sd_O2;                // [gr/day] O2
   
}

transformed parameters {

  // constraining the mean p_Am to be the same for all thermal grups
  array[N] real p_Am;

  {
    int start = 1;
    for (g in 1:n_temp) {
      int finish = start + n_group[g] - 1;
      real z_bar =
        mean(p_Am_star[start:finish]);
      real correction =
        sqrt(
          n_group[g] * 1.0 /
          (n_group[g] - 1.0)
        );
      for (i in start:finish) {
        p_Am[i] =
          p_Am_mu +
          p_Am_sd *
          correction *
          (p_Am_star[i] - z_bar);
      }
      start = finish + 1;
    }
  }
}

model {

  p_Am_mu ~ normal(prior_p_Am[1] , prior_p_Am[2]);
  p_Am_sd ~ normal(prior_p_Am_sd[1] , prior_p_Am_sd[2]);
  p_Am_star ~ normal(0.0,1.0);

  v ~ gamma(prior_v[1] , prior_v[2]);
  
  T_A ~ normal(prior_T_A[1] , prior_T_A[2]);
  T_A2 ~ normal(prior_T_A[1] , prior_T_A[2]);
  T_H ~ normal(prior_T_H[1] , prior_T_H[2]);
  T_AH ~ normal(prior_T_AH[1] , prior_T_AH[2]);
  T_L ~ normal(prior_T_L[1] , prior_T_L[2]);
  T_AL ~ normal(prior_T_AL[1] , prior_T_AL[2]);
  dH ~ normal (0.0,5.0);

  sd_length ~ normal(prior_sd_length[1],prior_sd_length[2]);
  sd_weight ~ normal(prior_sd_weight[1],prior_sd_weight[2]);
  sd_O2 ~ normal(prior_sd_O2[1],prior_sd_O2[2]);
    
  // Initial values for state variables
  array [N] vector [2] y0;
  for (i in 1:N){  // all fish (training and testing)
    E0_i[i] ~ normal(prior_E0_i[i], 2000.0);
    V0_i[i] ~ normal(prior_V0_i[i], 1.0);
    y0[i,1] = E0_i[i];
    y0[i,2] = V0_i[i];
    obs0[i,1] ~ normal( ((y0[i,2]^(1.0/3.0))/del_M), sd_length);              // length at t0
    obs0[i,2] ~ lognormal( log(y0[i,2] + y0[i,1]*w_E/(mu_E*d_V)), sd_weight); // weight at t0
  }

  // likelihood
  target += reduce_sum(partial_sum_lupdf,
                        obs[,,],   // array to be slice (multi threads_per_chain)
                        grainsize, // slicing grainsize (multi threads_per_chain)
                        // 1/0 filter for training fish
                        //is_train_fish[],
                        // arguments for ode_rk45 (DEB model)
                        y0[], t0[], age[,], n, N, 
                        // DEB main parameters
                        f,E_G,kap,p_Am[],p_M,k_J,v,kap_R,
                        // temperature-related paramters
                        Kelvin,Temp_mean[],Temp_amp[],pi2f,Temp_phi[],
                        T_A,T_A2,T_AL,T_L,T_AH,T_H,dH,T_ref,
                        birth_correction[],ID_temp[],Temperature_O2[],
                        // state variables to obervable variables 
                        del_M, d_V, w_E, mu_E, eta_O[,], n_M_n_O[,], 
                        // error
                        sd_length,sd_weight,sd_O2
                        );
  } 

generated quantities {
}

" # End of model code
,fill = TRUE)
sink()

#-----
# 3: Compiling model
#-----
#mod = cmdstan_model("DEB_O2.stan")
mod = cmdstan_model("DEB_O2.stan",cpp_options = list(stan_threads = TRUE)) # multicore settings

#-----
# 4: Initializing the model
#-----
# Setting priors (Supplementary materials 9)
V0_i=(obs0[,1]*parms$del_M)^(3)
prior_V0_i=V0_i
E0_i=(obs0[,2]-V0_i)*parms$mu_E*parms$d_V/parms$w_E
prior_E0_i=E0_i
prior_p_Am = c(parms$p_Am_after, 200)
prior_p_Am_sd = c(30,30)
prior_v = c(1.777778,8.888889)
prior_T_A=c(10000,5000) 
prior_T_A2=c(10000,5000)
prior_T_AH=c(100000,50000)
prior_T_H=c(305,5)
prior_T_AL=c(50000,20000)
prior_T_L=c(290,5)
prior_sd_length=c(0.3,0.1) 
prior_sd_weight=c(0.0,0.25)
prior_sd_O2=c(0.08,0.05)

n.chains = 4
initializer = function() list(
  "p_Am_mu"=prior_p_Am[1],
  "p_Am_sd"=prior_p_Am_sd[1],
  "p_Am_star"=rep(0,N),
  "v"=0.16,
  "E0_i"=prior_E0_i,
  "V0_i"=prior_V0_i,
  "sd_length"=prior_sd_length[1],
  "sd_weight"=0.01,
  "sd_O2"=prior_sd_O2[1],
  "T_A"=prior_T_A[1],
  "T_A2"=prior_T_A[1],
  "T_H"=prior_T_H[1],
  "T_AH"=prior_T_AH[1],
  "T_L"=prior_T_L[1],
  "T_AL"=prior_T_AL[1],
  "dH"=0.0
)
inits = list()
for (chain in 1:n.chains) inits[[chain]] = initializer()

#-----
# 5: Running
#-----

fit = mod$sample(
  data =list (
    # observations
    n = n,
    N = N,                                
    obs = obs_999,                            
    obs0 = obs0_999,
    t0 = t0,                          
    age = age,
    n_temp = n_temp,
    ID_temp =ID_temp,
    n_group=n_group,
    
    # fixed parameters
    f = parms$f,
    E_G = parms$E_G, 
    kap = parms$kap,
    p_M = parms$p_M,
    k_J = parms$k_J, 
    kap_R = parms$kap_R,
    del_M = parms$del_M,
    w_E = parms$w_E,        
    mu_E = parms$mu_E,
    d_V = parms$d_V,
    n_M_n_O = parms$n_M_n_O,
    eta_O = parms$eta_O,
    Kelvin = parms$Kelvin,
    Temp_mean = Temp_mean,
    Temp_amp = Temp_amp,
    pi2f = parms$pi2f,
    Temp_phi = Temp_phi,
    T_ref = parms$T_ref,
    birth_correction = birth_correction,
    Temperature_O2 = Temperature_O2,
    
    # Priors
    prior_E0_i=prior_E0_i,
    prior_V0_i=prior_V0_i,
    prior_p_Am=prior_p_Am,
    prior_p_Am_sd=prior_p_Am_sd,
    prior_v=prior_v,
    prior_sd_length=prior_sd_length,
    prior_sd_weight=prior_sd_weight,
    prior_sd_O2=prior_sd_O2,
    prior_T_A=prior_T_A,
    prior_T_H=prior_T_H,
    prior_T_AH=prior_T_AH,
    prior_T_L=prior_T_L,
    prior_T_AL=prior_T_AL,
    
    grainsize=4
  ),
  # settings
  chains = n.chains,
  parallel_chains = n.chains,
  threads_per_chain = 8,
  iter_warmup = 1000,
  iter_sampling = 5000,
  refresh = 100, 
  init = inits,
  adapt_delta = 0.80
)

#-----
# 6: Saving
#-----
fit$save_object(file = "out.RDS")
sink(file = "out_diagnostic.txt")
fit$cmdstan_diagnose()
sink()
