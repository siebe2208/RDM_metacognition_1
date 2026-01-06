library(tidyverse)
library(posterior)
library(bayesplot)
library(cmdstanr)


# -----------------------------
# Set parameters
# -----------------------------
n_trials <- 500        # number of trials
sigma_e <- 1    # noise on evidence (needs to be fixed for now)
sensitivity <-1.8 # perceptual sensitivity as a function of d' 

mu_beta <- 0.5               # criterion mean for the beta distribution
prec_beta <- 2       # criterion precision for the beta distribution
alpha <- mu_beta*prec_beta   # transformed shape parameter for beta distribution
beta <- (1-mu_beta)*prec_beta # transformed shape parameter for beta distribution

mu_k <- 0                  # criterion mean on the evidence scale (gaussian distribution)
sigma_k <- 1                 # criterion sd on the evidence scale (gaussian distribution)

rt_ndt = 0.2          # non decision time for RT
rt_sigma = 0.5        # Residual variance of RT
rt_e = 0.4            # scaling of the mean RT given evidence.

# Coherence level (stimulus intensity)
X <- abs(rnorm(n_trials,2,1))

# stimuli (left or right)
D <- sample(c(-1, 1), n_trials, replace = TRUE)  # true stimulus

# Simulate evidence
e <- rnorm(n_trials, mean = D*X, sd = sensitivity*sigma_e)

# Compute posterior P(D=1|e)
p_D1 = boot::inv.logit((2*e*1)/(sigma_e*sensitivity)^2)

# Posterior-based stochastic policy (using the beta distribution to set citeria)
c_beta <- extraDistr::rprop(n_trials,prec_beta,mu_beta)           # random threshold
p_beta <- pbeta(p_D1,alpha,beta) #probability of response on trial t given the citerion
a_beta <- ifelse(p_D1 > c_beta, 1, 0)         # choices

# Evidence-based stochastic policy (guassian distribution for citera)
k_gauss <- rnorm(n_trials, mu_k, sensitivity*sigma_k)          # citerion on trial t
p_gauss = pnorm((e-mu_k)/(sensitivity*sigma_k))                  #  probability of response on trial t given the citerion
a_gauss = rbinom(n_trials,1,p_gauss)              # binary response

# Dataframe
df <- data.frame(
  e = e,
  D = D,
  X =X,
  rt = rlnorm(n_trials,rt_e * -abs(e),rt_sigma) + rt_ndt,
  rt2 = rlnorm(n_trials,-abs(p_D1-0.5),rt_sigma) - rt_ndt,
  p_beta = p_beta,
  a_beta = a_beta,
  p_gauss = p_gauss,
  a_gauss = a_gauss)

# -----------------------------
# Plot simulations
# -----------------------------
df_long = df %>% pivot_longer(cols = c(a_beta, a_gauss, p_beta, p_gauss), 
                  names_to = c(".value", "noise_level"),  names_pattern = "(a|p)_(.*)") %>%
                  rename(resp = a, prob = p)

ggplot(df_long) + geom_point(aes(x=e, y=resp, color = noise_level), size=3)+ 
geom_line(aes(x=e, y=prob, color = noise_level), size =1.5)+
scale_x_continuous(limits = c(-max(abs(df_long$e)), max(abs(df_long$e))))+
theme_minimal()+
labs(x = "evidence", y = "P(choice=1)", color = "noise_level") +
scale_color_manual(values = c(beta = "salmon", gauss = "steelblue"))




