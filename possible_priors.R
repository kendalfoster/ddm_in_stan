
library("tidyverse")

n <- 1e5



### a
mu <- 0.75
s <- 0.5

df_a <- data.frame(
  t_dist = exp(extraDistr::rlst(n = n, df = 4, mu = 0.75, sigma = 0.75)),
  logis = exp(rlogis(n, location = mu, scale = s))
)

ggplot(data = df_a[df_a[["t_dist"]] < 100 & df_a[["logis"]] < 100, ]) +
  geom_histogram(aes(x = t_dist), binwidth = 0.25, boundary = 0, alpha = 0.7, fill = "blue") +
  geom_histogram(aes(x = logis), binwidth = 0.25, boundary = 0, alpha = 0.3, fill = "red") +
  geom_vline(xintercept = 1.25) +
  coord_cartesian(xlim = c(0, 10))



### ndt (t0)
mu <- -1
s <- .5

df_ndt <- data.frame(
  t_dist = exp(extraDistr::rlst(n = n, df = 4, mu = -1, sigma = 0.75)),
  logis = exp(rlogis(n, location = mu, scale = s))
)

ggplot(data = df_ndt[df_ndt[["t_dist"]] < 100 & df_ndt[["logis"]] < 100, ]) +
  geom_histogram(aes(x = t_dist), binwidth = 0.05, boundary = 0, alpha = 0.7, fill = "blue") +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.3, fill = "red") +
  geom_vline(xintercept = 0.2) +
  coord_cartesian(xlim = c(0, 1))



### w
mu <- 0
s <- 2/3

df_ndt <- data.frame(
  logis1 = plogis(rlogis(n = n, location = 0, scale = 0.65)),
  logis = plogis(rlogis(n, location = mu, scale = s))
)

ggplot(data = df_ndt[df_ndt[["logis1"]] < 100 & df_ndt[["logis"]] < 100, ]) +
  geom_histogram(aes(x = logis1), binwidth = 0.05, boundary = 0, alpha = 0.7, fill = "blue") +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.3, fill = "red") +
  geom_vline(xintercept = 0.5) +
  coord_cartesian(xlim = c(0, 1))



### sv
mu <- -0.5
s <- 0.5

df_ndt <- data.frame(
  t_dist = exp(extraDistr::rlst(n = n, df = 4, mu = -0.5, sigma = 0.7)),
  logis = exp(rlogis(n, location = mu, scale = s))
)

ggplot(data = df_ndt[df_ndt[["t_dist"]] < 100 & df_ndt[["logis"]] < 100, ]) +
  geom_histogram(aes(x = t_dist), binwidth = 0.05, boundary = 0, alpha = 0.7, fill = "blue") +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.3, fill = "red") +
  geom_vline(xintercept = 0.5) +
  coord_cartesian(xlim = c(0, 2.5))
