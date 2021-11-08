library("rstan")
model <- stan_model(file = file.path(getwd(), 'working_example.stan'),
                    verbose = TRUE)


##### check that the ddm_lpdf function can be used in fitting
library("rtdists")
N <- 100
a <- 1
v <- -2
t0 <- 0.3
# w <- 0.5 # rtdists uses z = w*a, defaults to a/2
sv <- 0.8
testd <- rtdists::rdiffusion(n = N, a = a, v = v, t0 = t0, sv = sv)
rt <- testd[["rt"]]
resp <- ifelse(testd[["response"]] == "upper", 1, 0)


m <- sampling(model, data = list(N=N, rt=rt, resp=resp))
m
