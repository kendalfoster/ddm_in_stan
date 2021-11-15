library("ggplot2")

n <- 1e5



### a
mu <- 0.75
s <- 0.5

ggplot(data = data.frame(logis = exp(rlogis(n, location = 0.75, scale = 0.5)))) +
  geom_histogram(aes(x = logis), binwidth = 0.25, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "a, threshold separation", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))


# df_a <- data.frame(
#   t_dist = exp(extraDistr::rlst(n = n, df = 4, mu = 0.75, sigma = 0.75)),
#   logis = exp(rlogis(n, location = mu, scale = s))
# )
# 
# ggplot(data = df_a[df_a[["t_dist"]] < 100 & df_a[["logis"]] < 100, ]) +
#   geom_histogram(aes(x = t_dist), binwidth = 0.25, boundary = 0, alpha = 0.7, fill = "blue") +
#   geom_histogram(aes(x = logis), binwidth = 0.25, boundary = 0, alpha = 0.3, fill = "red") +
#   geom_vline(xintercept = 1.25) +
#   coord_cartesian(xlim = c(0, 10))



### ndt (t0)
mu <- -1
s <- 0.5

ggplot(data = data.frame(logis = exp(rlogis(n, location = -1, scale = 0.5)))) +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "ndt, non-decision time", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# df_ndt <- data.frame(
#   t_dist = exp(extraDistr::rlst(n = n, df = 4, mu = -1, sigma = 0.75)),
#   logis = exp(rlogis(n, location = mu, scale = s))
# )
# 
# ggplot(data = df_ndt[df_ndt[["t_dist"]] < 100 & df_ndt[["logis"]] < 100, ]) +
#   geom_histogram(aes(x = t_dist), binwidth = 0.05, boundary = 0, alpha = 0.7, fill = "blue") +
#   geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.3, fill = "red") +
#   geom_vline(xintercept = 0.2) +
#   coord_cartesian(xlim = c(0, 1))



### w
mu <- 0
s <- 2/3

ggplot(data = data.frame(logis = plogis(rlogis(n, location = 0, scale = 2/3)))) +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "w, initial bias", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# df_ndt <- data.frame(
#   logis1 = plogis(rlogis(n = n, location = 0, scale = 0.65)),
#   logis = plogis(rlogis(n, location = mu, scale = s))
# )
# 
# ggplot(data = df_ndt[df_ndt[["logis1"]] < 100 & df_ndt[["logis"]] < 100, ]) +
#   geom_histogram(aes(x = logis1), binwidth = 0.05, boundary = 0, alpha = 0.7, fill = "blue") +
#   geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.3, fill = "red") +
#   geom_vline(xintercept = 0.5) +
#   coord_cartesian(xlim = c(0, 1))



### sv
mu <- -0.5
s <- 0.5

ggplot(data = data.frame(logis = exp(rlogis(n, location = -0.5, scale = 0.5)))) +
  geom_histogram(aes(x = logis), binwidth = 0.1, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 2.5)) +
  labs(x = "sv, variability in the drift rate", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

# df_ndt <- data.frame(
#   t_dist = exp(extraDistr::rlst(n = n, df = 4, mu = -0.5, sigma = 0.7)),
#   logis = exp(rlogis(n, location = mu, scale = s))
# )
# 
# ggplot(data = df_ndt[df_ndt[["t_dist"]] < 100 & df_ndt[["logis"]] < 100, ]) +
#   geom_histogram(aes(x = t_dist), binwidth = 0.05, boundary = 0, alpha = 0.7, fill = "blue") +
#   geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.3, fill = "red") +
#   geom_vline(xintercept = 0.5) +
#   coord_cartesian(xlim = c(0, 2.5))




# n <- 1e5
# 
# df_priors <- data.frame(
#   parameter = c(
#                 rep("a", n),
#                 rep("t0", n),
#                 rep("sv", n),
#                 rep("w", n)
#                 ),
#   value = c(
#             exp(rlogis(n, location = -1, scale = 0.5)),
#             exp(rlogis(n, location = 0.75, scale = 0.5)),
#             exp(rlogis(n, location = 0, scale = 2/3)),
#             exp(rlogis(n, location = -0.5, scale = 0.5))
#             )
#   )
# 
# ggplot(data = df_priors) +
#   geom_histogram(aes(x = value), binwidth = 0.1, boundary = 0, alpha = 0.7) +
#   coord_cartesian(xlim = c(0, 2.5)) +
#   labs(x = "sv, variability in the drift rate", y = "Count") +
#   theme_bw() +
#   theme(panel.border = element_blank(),
#         axis.text.x = element_text(size = 16),
#         axis.text.y = element_text(size = 16),
#         axis.title.x = element_text(size = 20),
#         axis.title.y = element_text(size = 20)) +
#   facet_wrap(~ factor(parameter, levels = c("a", "t0", "sv", "w")),
#              scales = "free_y", ncol = 2)



library("ggplot2")
library("cowplot")
n <- 1e5

p_a <- ggplot(data = data.frame(logis = exp(rlogis(n, location = 0.75, scale = 0.5)))) +
  geom_histogram(aes(x = logis), binwidth = 0.25, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 10)) +
  labs(x = "a, threshold separation", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

p_ndt <- ggplot(data = data.frame(logis = exp(rlogis(n, location = -1, scale = 0.5)))) +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "ndt, non-decision time", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

p_w <- ggplot(data = data.frame(logis = plogis(rlogis(n, location = 0, scale = 2/3)))) +
  geom_histogram(aes(x = logis), binwidth = 0.05, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = "w, initial bias", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

p_sv <- ggplot(data = data.frame(logis = exp(rlogis(n, location = -0.5, scale = 0.5)))) +
  geom_histogram(aes(x = logis), binwidth = 0.1, boundary = 0, alpha = 0.7) +
  coord_cartesian(xlim = c(0, 2.5)) +
  labs(x = "sv, variability in the drift rate", y = "Count") +
  theme_bw() +
  theme(panel.border = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))

cowplot::plot_grid(p_a, p_ndt, p_sv, p_w, ncol = 2)
