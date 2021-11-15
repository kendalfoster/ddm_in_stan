library("brms")
library("rtdists")
library("fddm")

ddm <- custom_family(
  "ddm", dpars = c("mu", "a", "ndt", "w", "sv"),
  links = c("identity", "log", "log", "logit", "log"), 
  lb = c(NA, 0, 0, 0, 0), ub = c(NA, NA, NA, 1, NA),
  type = "real", vars = "resp[n]"
)

stan_funs <- "
real ddm_lpdf(real rt, real v, real a, real ndt, real w, real sv, int resp)
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
  real v_i = (resp == 1) ? v : -v;
  real w_i = (resp == 1) ? w : 1 - w;
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
"

stanvars <- stanvar(scode = stan_funs, block = "functions")

log_lik_ddm <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  a <- brms::get_dpar(prep, "a", i = i)
  w <- brms::get_dpar(prep, "w", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  sv <- brms::get_dpar(prep, "sv", i = i)
  y <- prep$data$Y[i]
  resp <- prep$data$resp[i]
  fddm::dfddm(y, resp, a, mu, ndt, w, sv, log = TRUE)
}

posterior_predict_ddm <- function(i, prep, negative_rt = FALSE, ...) {
  if (!require("rtdists")) {
    stop("Package 'rtdists' must be installed for the function 'posterior_predict_ddm'.")
  }
  
  v <- as.numeric(brms::get_dpar(prep, "mu", i = i))
  a <- as.numeric(brms::get_dpar(prep, "a", i = i))
  w <- as.numeric(brms::get_dpar(prep, "w", i = i))
  t0 <- as.numeric(brms::get_dpar(prep, "ndt", i = i))
  sv <- as.numeric(brms::get_dpar(prep, "sv", i = i))
  n <- max(lengths(list(v, a, w, t0, sv)), 1) # I'm pretty sure it's always 1
  
  out <- rtdists::rdiffusion(
    n = n,
    a = a,
    v = v,
    t0 = t0,
    z = w * a,
    sv = sv
  )
  colnames(out) <- c("q", "resp")
  
  if (negative_rt) { # code "lower" responses as negative RTs
    out[["q"]] <- out[["q"]] * ifelse(out[["resp"]] == "upper", 1, -1)
  }
  
  return(out[["q"]])
}





data(speed_acc, package = "rtdists")
speed_acc <- droplevels(speed_acc[!speed_acc$censor,]) # remove extreme RTs
speed_acc <- droplevels(speed_acc[ speed_acc$frequency %in% 
                                     c("high", "nw_high"),])
speed_acc$resp <- as.integer(speed_acc$response)
speed_acc <- speed_acc[speed_acc[["id"]] == 1, ]


# library("bayestestR")
# options(contrasts = c('contr.orthonorm', 'contr.poly'))
# model.matrix(~condition*frequency, speed_acc)

# formula <- brmsformula(rt ~ 0 + condition:frequency +
#                             (0 + condition:frequency|p|id), 
#                        a ~ 0 + condition + (0 + condition|p|id), 
#                        ndt ~ 0 + condition + (0 + condition|p|id),
#                        w ~ 0 + condition + (0 + condition|p|id),
#                        sv ~ 1 + (1|id) )


get_prior(rt ~ 1,
          data = speed_acc,
          family = ddm)

# priors <- c(
#   set_prior("normal(0, 2.5)", class = "b"), # for mu (which is v)
#   set_prior("logistic(0.75, 0.5)", class = "b", dpar = "a"),
#   set_prior("logistic(-1, 0.5)", class = "b", dpar = "ndt"),
#   set_prior("logistic(0, 0.67)", class = "b", dpar = "w"),
#   set_prior("logistic(-0.5, 0.5)", class = "sd", dpar = "sv")
# )
priors <- c(
  set_prior("normal(0, 2.5)", class = "Intercept"), # for mu (which is v)
  set_prior("logistic(0.75, 0.5)", class = "a"),
  set_prior("logistic(-1, 0.5)", class = "ndt"),
  set_prior("logistic(0, 0.67)", class = "w"),
  set_prior("logistic(-0.5, 0.5)", class = "sv")
)

# make_stancode(formula,
#               family = ddm,
#               data = speed_acc,
#               prior = priors)
# 
# 
# tmp_dat <- make_standata(formula,
#                          family = ddm,
#                          data = speed_acc,
#                          prior = priors)
# str(tmp_dat)
# init_fun <- function() {
#   list(
#     b = rnorm(tmp_dat$K),
#     b_a = exp(runif(tmp_dat$K_a, 1, 2)),
#     b_ndt = exp(runif(tmp_dat$K_ndt, 0.1, 0.15)),
#     b_w = exp(rnorm(tmp_dat$K_w, 0.5, 1)),
#     b_sv = exp(runif(tmp_dat$K_sv, 0.1, 0.7))
#     # sd_1 = runif(tmp_dat$M_1, 0.5, 1),
#     # z_1 = matrix(rnorm(tmp_dat$M_1*tmp_dat$N_1, 0, 0.01),
#     #              tmp_dat$M_1, tmp_dat$N_1),
#     # L_1 = diag(tmp_dat$M_1),
#     # sd_2 = runif(tmp_dat$M_2, 0.5, 1),
#     # z_2 = matrix(rnorm(tmp_dat$M_2*tmp_dat$N_2, 0, 0.01),
#     #              tmp_dat$M_2, tmp_dat$N_2)
#   )
# }


stanvars_rw <- stanvars +
  stanvar(x = as.integer(speed_acc[["resp"]]), 
          name = "resp", scode = "int resp[N];")

fit_rwdata_5par <- brm(rt ~ 1,
                       family = ddm,
                       prior = priors,
                       # inits = init_fun,
                       stanvars = stanvars_rw,
                       data = speed_acc,
                       chains = 1,
                       # iter = 500,#2000
                       # warmup = 250,#1000
                       thin = 1,
                       cores = getOption("mc.cores", 1),
                       control = list(
                         max_treedepth = 10,
                         adapt_delta = 0.8
                       )
)

# saveRDS(fit_rwdata_5par,
#         file = file.path(getwd(),"fits", "fit_rwdata_5par.RDS"))
# 
# fit_rwdata_5par <- readRDS(file = file.path(getwd(), "fits",
#                                             "fit_rwdata_5par.RDS"))

