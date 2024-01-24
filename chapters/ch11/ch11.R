library(rethinking)
library(tidyverse)

# Confound generative simulation ------------------------------------------

set.seed(12)
N <- 2000
G <- sample(1:2, N, replace = TRUE)
# "ability" high (1) or average (0)
u <- rbern(N, 0.1)
# gender 1 tends to apply to dept 1, and 2 to 2
# G=1 w/ higher ability tends to apply to 2
D <- rbern(N, ifelse(G == 1, u * 1, 0.75)) + 1

table(G, D, u)

# matrix of acceptance rates
p_u0 <- matrix(c(0.1, 0.1, 0.1, 0.3), nrow = 2)
p_u1 <- matrix(c(0.3, 0.3, 0.5, 0.5), nrow = 2)
p_u <- list(p_u0, p_u1)
# simulate acceptance
p <- sapply(seq(N), \(i) p_u[[1 + u[i]]][D[i], G[i]])
A <- rbern(N, p)

dat_sim <- list(A = A, G = G, D = D)
# total effect of gender
m1 <- ulam(
  alist(
    A ~ bernoulli(p),
    logit(p) <- a[G],
    a[G] ~ normal(0, 1)
  ), data = dat_sim, chains = 4, cores = 4
)
# direct effects - now confounded!
m2 <- ulam(
  alist(
    A ~ bernoulli(p),
    logit(p) <- a[G, D],
    matrix[G, D]:a ~ normal(0, 1)
  ), data = dat_sim, chains = 4, cores = 4
)

# total effect of gender is negative
precis(m1, depth = 3)
# indirect effect of gender is confounded
precis(m2, depth = 3)
post2 <- extract.samples(m2)
# D1G1
dens(inv_logit(post2$a[, 1, 1]), 
     lwd = 3, col = 4, xlim = c(0, 0.5), xlab = "probability of admission")
# D1G2
dens(inv_logit(post2$a[, 2, 1]), 
     lwd = 3, col = 4, lty = 2, add = TRUE)
# D2G1
dens(inv_logit(post2$a[, 1, 2]), 
     lwd = 3, col = 2, add = TRUE)
# D2G2
dens(inv_logit(post2$a[, 2, 2]), 
     lwd = 3, col = 2, lty = 2, add = TRUE)


# Sensitivity analysis ----------------------------------------------------

datl <- dat_sim
datl$D2 <- ifelse(datl$D == 2, 1, 0)
datl$N <- length(datl$D)
datl$b <- c(1, 1)
datl$g <- c(1, 0)

mGDu <- ulam(
  alist(
    # A model
    A ~ bernoulli(p),
    logit(p) <- a[G, D] + b[G] * u[i],
    matrix[G, D]:a ~ normal(0, 1),
    
    # D2 model
    D2 ~ bernoulli(q),
    logit(q) <- delta[G] + g[G] * u[i],
    delta[G] ~ normal(0, 1),
    
    # declare unobserved u
    vector[N]:u ~ normal(0, 1)
  ), data = datl, chains = 4, cores = 4
)

postGDu <- extract.samples(mGDu)
# D1G1
dens(inv_logit(postGDu$a[, 1, 1]), 
     lwd = 3, col = 4, xlim = c(0, 0.5), xlab = "probability of admission")
# D1G2
dens(inv_logit(postGDu$a[, 2, 1]), 
     lwd = 3, col = 4, lty = 2, add = TRUE)
# D2G1
dens(inv_logit(postGDu$a[, 1, 2]), 
     lwd = 3, col = 2, add = TRUE)
# D2G2
dens(inv_logit(postGDu$a[, 2, 2]), 
     lwd = 3, col = 2, lty = 2, add = TRUE)

# Oceania tech ------------------------------------------------------------

data(Kline)
d <- Kline
d$P <- scale(log(d$population))
d$contact_id <- ifelse(d$contact == "high", 2, 1)
dat <- list(
  T = d$total_tools,
  P = d$P,
  C = d$contact_id
)
# intercept only
m11.9 <- ulam(
  alist(
    T ~ dpois(lambda),
    log(lambda) <- a,
    a ~ dnorm(3, 0.5)
  ), 
  data = dat, chains = 4, log_lik = TRUE
)
# interaction model
m11.10 <- ulam(
  alist(
    T ~ dpois(lambda),
    log(lambda) <- a[C] + b[C] * P,
    a[C] ~ dnorm(3, 0.5),
    b[C] ~ dnorm(0, 0.2)
  ), 
  data = dat, chains = 4, log_lik = TRUE
)
compare(m11.9, m11.10, func = PSIS)

# Visualize
k <- PSIS(m11.10, pointwise = TRUE)$k
plot(dat$P, dat$T, 
     xlab = "log population (st)", ylab = "total tools",
     col = ifelse(dat$C == 1, 4, 2),
     lwd = 4 + 4 * normalize(k),
     ylim = c(0, 75),
     cex = 1 + normalize(k))
# set up x-axis values for predictions
P_seq <- seq(from = -1.4, to = 3, len = 100)

# Predictions for C=1 (low contact)
lambda <- link(m11.10, data = data.frame(P = P_seq, C = 1))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty = 2, lwd = 1.5)
shade(lci, P_seq, xpd = TRUE, col = col.alpha(4, 0.3))

# Predictions for C=2 (high contact)
lambda <- link(m11.10, data = data.frame(P = P_seq, C = 2))
lmu <- apply(lambda, 2, mean)
lci <- apply(lambda, 2, PI)
lines(P_seq, lmu, lty = 1, lwd = 1.5)
shade(lci, P_seq, xpd = TRUE, col = col.alpha(2, 0.3))

# Tool equilibria
f <- function(a = 0.02, b = 0.5, g = 0.2, P = 1e4, t_max = 50) {
  T <- rep(0, t_max)
  for (i in 2:t_max)
    T[i] <- T[i - 1] + a * P^b - g * T[i - 1]
  return (T)
}

plot(NULL, 
     xlim = c(0, 50), ylim = c(0, 10),
     xlab = "time", ylab = "tools")
T <- f(P = 1e3)
lines(1:50, T, lwd = 3, col = 2)
T <- f(P = 1e4)
lines(1:50, T, lwd = 3, col = 2)

# innovation/loss model

dat2 <- list(
  T = d$total_tools,
  P = d$population,
  C = d$contact_id
)
m11.11 <- ulam(
  alist(
    T ~ dpois(lambda),
    lambda <- exp(a[C]) * P^b[C] / g,
    a[C] ~ dnorm(1, 1),
    b[C] ~ dexp(1),
    g ~ dexp(1)
  ),
  data = dat2, chains = 4, cores = 4
)
