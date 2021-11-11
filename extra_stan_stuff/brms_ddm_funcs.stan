// vector inputs
real ddm_lpdf(vector rt, real v, real a, real ndt, real w, real sv, vector resp)
{
  // initialize output
  real logp = 0.0;
  
  // get lengths of inputs
  int nrt = num_elements(rt);
  int nresp = num_elements(resp);
  vector[nrt] new_resp;
  if (nrt != nresp) {
    for (i in 1:nrt) {
      new_resp[i] = resp[1];
    }
  } else {
    for (i in 1:nresp) {
      new_resp[i] = resp[i];
    }
  }
  
  // check parameter values
  if (a <= 0 || is_inf(a) || is_nan(a)) return negative_infinity();
  if (is_inf(v) || is_nan(v)) return negative_infinity();
  if (ndt < 0 || is_inf(ndt) || is_nan(ndt)) return negative_infinity();
  if (w <= 0 || w >= 1 || is_nan(w)) return negative_infinity();
  if (sv < 0 || is_inf(sv) || is_nan(sv)) return negative_infinity();
  
  // loop through inputs
  for (i in 1:nrt) {
    // adjust parameters
    real t = rt[i] - ndt;
    real a_i = a;
    real sv_i = sv;
    real v_i = (new_resp[i] < 0.5) ? v : -v;
    real w_i = (new_resp[i] < 0.5) ? w : 1 - w;
    real taa = t / (a_i * a_i);
    
    // define variables for use later (have to do it early because of Stan)
    // constants
    real log_2_pi_2 = 0.5 * log(2 * pi());
    real sv_thresh = 0.0;
    // used in both large-time and small-time
    real mult = 0.0;
    real exp_err = 0.0;
    real gamma = -0.5 * pi() * pi() * taa;
    real summ = 0.0;
    int max_int = 1000;
    // used only in large-time
    int kl = 0;
    real k_dbl = 0.0;
    real bc = 1 / (pi() * sqrt(taa));
    // used only in small-time
    real minterms = 0.5 * sqrt(taa) - 0.5 * w_i; // min number of terms, truncates toward 0
    real term = 0.0;
    real rj = 0.0;
    int js = 0;
    
    // check large-time number of terms
    if (sv_i < sv_thresh) {
      mult = - v_i * a_i * w_i - 0.5 * v_i*v_i * t - 2 * log(a_i);
    } else {
      mult = (sv_i*sv_i * a_i*a_i * w_i*w_i - 2 * v_i * a_i * w_i - v_i*v_i * t)
      / (2 + 2 * sv_i*sv_i * t) - 0.5 * log(1 + sv_i*sv_i * t) - 2 * log(a_i);
    }
    exp_err = 0.000001 * exp(-mult);
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
      if (sv_i < sv_thresh) {
        mult = log(a_i) - log_2_pi_2 - 1.5 * log(t) - v_i * a_i * w_i
        - 0.5 * v_i*v_i * t;
      } else {
        mult = log(a_i) - 1.5 * log(t) - log_2_pi_2
        - 0.5 * log(1 + sv_i*sv_i * t)
        + (sv_i*sv_i * a_i*a_i * w_i*w_i - 2 * v_i * a_i * w_i
        - v_i*v_i * t) / (2 + 2 * sv_i*sv_i * t);
      }
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
      rj = 2 * js - w_i;
      term = rj * exp(gamma * rj*rj);
      summ -= term;
      while (term > exp_err) {
        rj = 2 * js + w_i;
        term = rj * exp(gamma * rj*rj);
        summ += term;
        if (term <= exp_err) break;
        js += 1;
        rj = 2 * js - w_i;
        term = rj * exp(gamma * rj*rj);
        summ -= term;
      }
    }
    
    // check summ and chuck everything in logp
    if (summ >= 0) {  // if result is negative, don't add to logp
      logp += mult + log(summ);
    }
  }
  
  return logp;
}

// non-vector (i.e., real) inputs
real ddm_lpdf(real rt, real v, real a, real ndt, real w, real sv, real resp)
{
  // initialize output
  real logp = 0.0;
  
  // check parameter values
  if (a <= 0 || is_inf(a) || is_nan(a)) return negative_infinity();
  if (is_inf(v) || is_nan(v)) return negative_infinity();
  if (ndt < 0 || is_inf(ndt) || is_nan(ndt)) return negative_infinity();
  if (w <= 0 || w >= 1 || is_nan(w)) return negative_infinity();
  if (sv < 0 || is_inf(sv) || is_nan(sv)) return negative_infinity();
  
  // adjust parameters
  {
  real t = rt - ndt;
  real a_i = a;
  real sv_i = sv;
  real v_i = (resp < 0.5) ? v : -v;
  real w_i = (resp < 0.5) ? w : 1 - w;
  real taa = t / (a_i * a_i);
  
  // define variables for use later (have to do it early because of Stan)
  // constants
  real log_2_pi_2 = 0.5 * log(2 * pi());
  real sv_thresh = 0.0;
  // used in both large-time and small-time
  real mult = 0.0;
  real exp_err = 0.0;
  real gamma = -0.5 * pi() * pi() * taa;
  real summ = 0.0;
  int max_int = 1000;
  // used only in large-time
  int kl = 0;
  real k_dbl = 0.0;
  real bc = 1 / (pi() * sqrt(taa));
  // used only in small-time
  real minterms = 0.5 * sqrt(taa) - 0.5 * w_i; // min number of terms, truncates toward 0
  real term = 0.0;
  real rj = 0.0;
  int js = 0;
  
  // check large-time number of terms
  if (sv_i < sv_thresh) {
    mult = - v_i * a_i * w_i - 0.5 * v_i*v_i * t - 2 * log(a_i);
  } else {
    mult = (sv_i*sv_i * a_i*a_i * w_i*w_i - 2 * v_i * a_i * w_i - v_i*v_i * t)
    / (2 + 2 * sv_i*sv_i * t) - 0.5 * log(1 + sv_i*sv_i * t) - 2 * log(a_i);
  }
  exp_err = 0.000001 * exp(-mult);
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
    if (sv_i < sv_thresh) {
      mult = log(a_i) - log_2_pi_2 - 1.5 * log(t) - v_i * a_i * w_i
      - 0.5 * v_i*v_i * t;
    } else {
      mult = log(a_i) - 1.5 * log(t) - log_2_pi_2
      - 0.5 * log(1 + sv_i*sv_i * t)
      + (sv_i*sv_i * a_i*a_i * w_i*w_i - 2 * v_i * a_i * w_i
      - v_i*v_i * t) / (2 + 2 * sv_i*sv_i * t);
    }
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
    rj = 2 * js - w_i;
    term = rj * exp(gamma * rj*rj);
    summ -= term;
    while (term > exp_err) {
      rj = 2 * js + w_i;
      term = rj * exp(gamma * rj*rj);
      summ += term;
      if (term <= exp_err) break;
      js += 1;
      rj = 2 * js - w_i;
      term = rj * exp(gamma * rj*rj);
      summ -= term;
    }
    
    // check summ and chuck everything in logp
    if (summ >= 0) {  // if result is negative, don't add to logp
      logp += mult + log(summ);
    }
  }
  }
  
  return logp;
}
