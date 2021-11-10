library("rstan")
model <- stan_model(file = file.path(getwd(), 'working_example.stan'),
                    verbose = FALSE)


##### check that the ddm_lpdf function can be used in fitting
library("rtdists")
N <- 1000
a <- 1
v <- -2
t0 <- 0.3
# w <- 0.5 # rtdists uses z = w*a, defaults to a/2
sv <- 0.8
set.seed(1234567)
testd <- rtdists::rdiffusion(n = N, a = a, v = v, t0 = t0, sv = sv)
rt <- testd[["rt"]]
resp <- ifelse(testd[["response"]] == "upper", 1, 0)


m <- sampling(model,
              chains = 1,
              data = list(N = N,
                          rt = rt,
                          resp = resp),
              init = list(list(
                               a = 1.1,
                               v = -0.2,
                               t0 = 0.2,
                               w = 0.45,
                               sv = 0.7
                               )))
m




##### check with MLE
library("fddm")

ll_fun <- function(pars, rt, resp, err_tol) {
  dens <- dfddm(rt = rt, response = resp, a = pars[[1]], v = pars[[2]], t0 = pars[[3]],
                w = pars[[4]], sv = pars[[5]], log = TRUE, err_tol = 1e-6)
  
  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

fit <- nlminb(c(1.1, -0.2, 0.2, 0.45, 0.05), objective = ll_fun,
              rt = rt, resp = resp + 1,
              # limits:   a,    v,  t0, w,  sv
              lower = c(.01, -Inf,   0, 0,   0),
              upper = c(Inf,  Inf, Inf, 1, Inf))
fit

