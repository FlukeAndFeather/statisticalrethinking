library(rethinking)
library(tidyverse)
data(UCBadmit)
d <- UCBadmit

# Note: corresponds to lecture 09

# Berkeley admissions: data exploration -----------------------------------

# generative model, basic mediator scenario
N <- 1000
# Gender
G <- sample(1:2, size = N, replace = TRUE)
# Gender 1 tends to apply to department, 2 to 2
D <- rbern(N, ifelse(G == 1, 0.3, 0.8)) + 1
# matrix of acceptance rates [dept, gender]
accept_rate <- matrix(c(0.1, 0.3, 0.1, 0.3), nrow = 2)
# simulate acceptance
A <- rbern(N, accept_rate[D, G])

table(G, D) # more Gender 1 -> Dept 1, Gender 2 -> Dept 2
table(G, A)
# accept rates ~0.19, 0.27, depending on simulation
tibble(G, A) %>% 
  group_by(G) %>% 
  summarize(accept_rate = mean(A))

# Introduce direct discrimination
accept_rate <- matrix(c(0.05, 0.2, 0.1, 0.3), nrow = 2)
A <- rbern(N, accept_rate[D, G])

table(G, D) # more Gender 1 -> Dept 1, Gender 2 -> Dept 2
table(G, A)
# accept rates ~0.09, 0.16, depending on simulation
tibble(G, A) %>% 
  group_by(G) %>% 
  summarize(accept_rate = mean(A))

# Overall, same pattern as absence of direct discrimination!

# Logistic modeling -------------------------------------------------------

# "Flat" prior
a <- rnorm(1e4, 0, 10)
b <- rnorm(1e4, 0, 10)

xseq <- seq(from = -3, to = 3, length.out = 100)
p <- sapply(xseq, function(x) inv_logit(a + b * x))

plot(NULL, 
     xlim = c(-2.5, 2.5), ylim = c(0, 1),
     xlab = "x value", ylab = "probability")
for (i in 1:10) {
  lines(xseq, p[i, ], lwd = 3, col = 2)
}

# look how strong those relationships are!

# Compare to narrow prior. Assumes less of a relationship
a <- rnorm(1e4, 0, 1.5)
b <- rnorm(1e4, 0, 0.5)

xseq <- seq(from = -3, to = 3, length.out = 100)
p <- sapply(xseq, function(x) inv_logit(a + b * x))

plot(NULL, 
     xlim = c(-2.5, 2.5), ylim = c(0, 1),
     xlab = "x value", ylab = "probability")
for (i in 1:10) {
  lines(xseq, p[i, ], lwd = 3, col = 2)
}

# Berkeley admissions: statistical model ----------------------------------

# generative model, basic mediator scenario
N <- 1000
# Gender
G <- sample(1:2, size = N, replace = TRUE)
# Gender 1 tends to apply to department, 2 to 2
D <- rbern(N, ifelse(G == 1, 0.3, 0.8)) + 1
# matrix of acceptance rates [dept, gender]
accept_rate <- matrix(c(0.05, 0.2, 0.1, 0.3), nrow = 2)
# simulate acceptance
A <- rbern(N, accept_rate[D, G])

dat_sim <- list(A = A, D = D, G = G)

m1 <- ulam(
  alist(
    A ~ bernoulli(p),
    logit(p) <- a[G],
    a[G] ~ normal(0, 1)
  ),
  data = dat_sim,
  chains = 4, 
  cores = 4
)

m2 <- ulam(
  alist(
    A ~ bernoulli(p),
    logit(p) <- a[G, D],
    matrix[G, D]:a ~ normal(0, 1)
  ),
  data = dat_sim,
  chains = 4, 
  cores = 4
)

# On the log-odds scale. Note param estimates <0 i.e. prob < 50%.
precis(m1, depth = 2)
# m2 conditions acceptance on department, too.
precis(m2, depth = 3)
# Coefs as probability. a[1:2,1] ~10%, a[1:2,2] ~30%
inv_logit(coef(m2))

m2_post <- extract.samples(m2)
a_post <- expand_grid(i = 1:2, j = 1:2) %>% 
  mutate(ij = paste(i, j, sep = ","),
         a = map2(i, j, \(i, j) m2_post$a[, i, j])) %>% 
  unnest(a) %>% 
  mutate(a_prob = inv_logit(a))
a_post %>% 
  group_by(ij) %>% 
  summarize(mean_a = mean(a)) %>% 
  pull(mean_a) %>% 
  matrix(nrow = 2)
ggplot(a_post, aes(a_prob, color = ij)) +
  geom_density() +
  theme_classic()

# Berkeley admissions: analyze --------------------------------------------

dat <- list(
  A = UCBadmit$admit,
  N = UCBadmit$applications,
  G = ifelse(UCBadmit$applicant.gender == "female", 1, 2),
  D = as.integer(UCBadmit$dept)
)

# total effect of gender
mG <- ulam(
  alist(
    A ~ binomial(N, p),
    logit(p) <- a[G],
    a[G] ~ normal(0, 1)
  ),
  data = dat, chains = 4, cores = 4
)

mGD <- ulam(
  alist(
    A ~ binomial(N, p),
    logit(p) <- a[G, D],
    matrix[G, D]:a ~ normal(0, 1)
  ),
  data = dat, chains = 4, cores = 4
)

traceplot(mG)
traceplot(mGD)
 
precis(mG, depth = 2)
precis(mGD, depth = 3)

mG_post <- extract.samples(mG)
prA_G1 <- inv_logit(mG_post$a[, 1])
prA_G2 <- inv_logit(mG_post$a[, 2])
mG_contrast <- prA_G1 - prA_G2
dens(mG_contrast, lwd = 4, col = 2, xlab = "Gender contrast (probability)")

mGD_post <- extract.samples(mGD)
prA <- inv_logit(mGD_post$a)
mGD_a_contrast <- sapply(1:6, \(i) prA[, 1, i] - prA[, 2, i])
plot(NULL, 
     xlim = c(-0.2, 0.3), ylim = c(0, 25), 
     xlab = "Gender contrast (probability)",
     ylab = "Density")
for (i in 1:6) {
  dens(mGD_a_contrast[, i], lwd = 4, col = 1 + i, add = TRUE)
}

# Berkeley admissions: intervention ---------------------------------------

# simulate applications per department
total_apps <- sum(dat$N)
apps_per_dept <- sapply(1:6, \(i) sum(dat$N[dat$D == i]))

# simulate as if all apps men or women
p_G1 <- link(mGD, 
             data = list(
               D = rep(1:6, times = apps_per_dept),
               N = rep(1, total_apps),
               G = rep(1, total_apps)
             ))
p_G2 <- link(mGD, 
             data = list(
               D = rep(1:6, times = apps_per_dept),
               N = rep(1, total_apps),
               G = rep(2, total_apps)
             ))

# summarize
dev.off()
dens(p_G1 - p_G2, lwd = 4, col = 2, xlab = "effect of gender perception")
abline(v = 0, lty = 2)
text(-0.05, 8, "men adv", adj = c(1, 1), font = 3)
text(0.05, 8, "women adv", adj = c(0, 1), font = 3)

table(p_G1 - p_G2 > 0)
