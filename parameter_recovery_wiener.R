library("brms")

par_rec_dataw <- readRDS(file.path(getwd(), "fits",
                                  "parameter_recovery_data.RDS"))
par_rec_dataw[["response2"]] <- as.integer(par_rec_dataw[["response"]]) - 1

get_prior(rt | dec(response2) ~ 1,
          family = wiener(link_bs = "log",
                          link_ndt = "log",
                          link_bias = "logit"),
          data = par_rec_data)

priorsw <- c(
  set_prior("normal(0, 2.5)", class = "Intercept"), # for mu (which is v)
  set_prior("logistic(0.75, 0.5)", class = "bs"),
  set_prior("logistic(-1, 0.5)", class = "ndt"),
  set_prior("logistic(0, 0.65)", class = "bias")
)

fit_parrec_wiener <- brm(rt | dec(response2) ~ 1,
                         family = wiener(link_bs = "log",
                                         link_ndt = "log",
                                         link_bias = "logit"),
                         prior = priorsw,
                         data = par_rec_dataw,
                         chains = 1,
                         iter = 2000,
                         warmup = 1000,
                         thin = 1,
                         cores = getOption("mc.cores", 1),
                         control = list(
                           max_treedepth = 10,
                           adapt_delta = 0.8
                           )
                         )
saveRDS(fit_parrec_wiener,
        file = file.path(getwd(), "fits", "fit_parrec_wiener.RDS"))

fit_parrec_wiener <- readRDS(file = file.path(getwd(), "fits",
                                              "fit_parrec_wiener.RDS"))