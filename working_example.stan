functions{
  real ddm_lpdf(real rt, real resp, real a, real v, real t0, real w, real sv)
  {
    // adjust parameters
    real t = rt - t0;
    real a_i = a;
    real sv_i = sv;
    real v_i = (resp < 0.5) ? v : -v;
    real w_i = (resp < 0.5) ? w : 1 - w;
    real taa = t / (a_i * a_i);
    
    // define variables for use later (have to do it early because of Stan)
    real log_2_pi_2 = 0.5 * log(2 * pi());
    real gamma = -0.5 * pi() * pi() * taa;
    real summ = 0.0;
    real minterms = 0.5 * sqrt(taa) - 0.5 * w_i; // min number of terms, truncates toward 0
    real term = 0.0;
    real rj = 0.0;
    int js = 0;

    // initialize output
    real logp = 0.0;
    
    // check large-time number of terms
    real mult = (sv_i * sv_i * a_i * a_i * w_i * w_i - 2 * v_i * a_i * w_i - v_i * v_i * t)
                / (2 + 2 * sv_i * sv_i * t) - 0.5 * log(1 + sv_i * sv_i * t) - 2 * log(a_i);
    real exp_err = 0.000001 * exp(-mult);
    int kl = 0;
    int max_int = 1000;
    real k_dbl = 0.0;
    real bc = 1 / (pi() * sqrt(taa));
    if (bc > max_int) { // boundary condition
      k_dbl = max_int;
    } else {
      if (exp_err * pi() * taa < 1) { // error threshold is low enough
        k_dbl = sqrt(-2 * log(pi() * taa * exp_err) / (pi() * pi() * taa));
        if (k_dbl > max_int) {
          k_dbl = max_int;
        } else {
          k_dbl = fmax(k_dbl, bc);
        }
      } else { // threshold not low enough, so set to boundary condition
        k_dbl = bc;
      }
    }
    while (kl < k_dbl) kl += 1;
    
    // compare large-time and small-time
    if (kl <= 1) { // do large-time
      for(jl in 1:kl) {
        summ += jl * sin(jl * w_i * pi()) * exp(gamma * jl * jl);
      }
      logp += log(pi());
    } else { // do small-time
      mult = log(a_i) - 1.5 * log(t) - log_2_pi_2 - 0.5 * log(1 + sv_i*sv_i * t)
             + (sv_i*sv_i * a_i*a_i * w_i*w_i - 2 * v_i * a_i * w_i
                - v_i*v_i * t) / (2 + 2 * sv_i*sv_i * t);
      exp_err = 0.000001 * exp(-mult);
      gamma = -1 / (2 * taa);
      summ = w_i * exp(gamma * w_i*w_i); // initialize with j=0 term
      while (js < minterms) {
        js += 1;
        rj = 2 * js - w_i;
        summ -= rj * exp(gamma * rj*rj);
        rj = 2 * js + w_i;
        summ += rj * exp(gamma * rj*rj);
      }
      js += 1;
      rj = 2 * js - w;
      term = rj * exp(gamma * rj*rj);
      summ -= term;
      while (term > exp_err) {
        rj = 2 * js + w;
        term = rj * exp(gamma * rj*rj);
        summ += term;
        if (term <= exp_err) break;
        js += 1;
        rj = 2 * js - w;
        term = rj * exp(gamma * rj*rj);
        summ -= term;
      }
    }
    
    // check summ and chuck everything in logp
    if (summ >= 0) {  // if result is negative, don't add to logp
      logp += mult + log(summ);
    }
    
    return logp;
  }
}

data {
  int<lower=0> N; // number of observations
  vector[N] rt; // response time
  vector[N] resp; // response
}

parameters {
  real<lower=0> a;
  real v;
  real<lower=0> t0;
  real<lower=0, upper=1> w;
  real<lower=0> sv;
}

model {
  for(i in 1:N)
    rt[i] ~ ddm(resp[i], a, v, t0, w, sv);
}
