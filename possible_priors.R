
library("tidyverse")

n <- 1e5


### a
tib <- tibble(
  sample = seq_len(n),
  #linear = rnorm(n, 0.75, 0.7)
  linear = extraDistr::rlst(n = n, df = 4, 
                            mu = 0.75, sigma = 0.75)
) %>% 
  mutate(
    a = exp(linear)
  )

tib %>% 
  filter(a < 100) %>% 
  ggplot(aes(a)) +
  geom_histogram(binwidth = 0.25, boundary = 0) +
  geom_vline(xintercept = 1) +
  coord_cartesian(xlim = c(0, 10))

### ndt
tib <- tibble(
  sample = seq_len(n),
  #linear = rnorm(n, 0.75, 0.7)
  linear = extraDistr::rlst(n = n, d = 4, 
                            mu = -1, sigma = 0.75)
) %>% 
  mutate(
    b = exp(linear)
  )

tib %>% 
  filter(b < 100) %>% 
  ggplot(aes(b)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  geom_vline(xintercept = 0.2) +
  coord_cartesian(xlim = c(0, 1))

### w
tib <- tibble(
  sample = seq_len(n),
  #linear = rnorm(n, 0.75, 0.7)
  linear = rlogis(n = n, location = 0, 0.65)
) %>% 
  mutate(
    a = plogis(linear)
  )

tib %>% 
  filter(a < 100) %>% 
  ggplot(aes(a)) +
  geom_histogram(binwidth = 0.01, boundary = 0) +
  geom_vline(xintercept = 0.5) +
  coord_cartesian(xlim = c(0, 1))

### sv

tib <- tibble(
  sample = seq_len(n),
  #linear = rnorm(n, 0.75, 0.7)
  linear = extraDistr::rlst(n = n, d = 4, 
                            mu = -0.5, sigma = 1/sqrt(2))
) %>% 
  mutate(
    b = exp(linear)
  )

tib %>% 
  filter(b < 100) %>% 
  ggplot(aes(b)) +
  geom_histogram(binwidth = 0.05, boundary = 0) +
  geom_vline(xintercept = 0.5) +
  coord_cartesian(xlim = c(0, 2.5))

