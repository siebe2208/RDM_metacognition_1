library(cmdstanr)
library(here)
library(ggplot2)

data = read.csv(choose.files())
variances <- tapply(data$coherence, gl(5, 10), var)
print(variances)
hist(data$coherence[11:50])

main_data = data[data$Trialtype == "Main",]

main_data$evidence = ifelse(main_data$dots.direction == 180, main_data$coherence, -(main_data$coherence))
main_data$resp_right = ifelse(main_data$resp == "['left']", 0, 1)

x = main_data$evidence
y = main_data$resp_right

data_list = list(N = length(x), x=x, y=y)

model = cmdstan_model("model_RDM_1.stan")

fit = model$sample(data = data_list)

draws = fit$draws(c("b", "lapse"))

data_a = as.vector(fit$draws("a"))
hist(data_a)
mean_a = mean(data_a)

data_b = as.vector(fit$draws("b"))
hist(data_b)

mean_b = mean(data_b)
print(mean_b)

data_lapse = as.vector(fit$draws("lapse"))
hist(data_lapse)
mean_lapse=mean(data_lapse)

psycho = function(a,b,lapse,x){lapse + (1-2*lapse)*(1/(1+exp(-b*(x-a))))}

x_seq = seq(-1, 1, length.out = 200)

curve = psycho(mean_a, mean_b, mean_lapse, x_seq)

df_curve <- data.frame(x = x_seq, theta = curve)
df_data <- data.frame(x = x, y = y)

ggplot() +
  geom_line(data = df_curve, aes(x = x, y = theta), color = "blue", size = 1.2) +
  geom_point(data = df_data, aes(x = x, y = y), color = "red") +
  ylim(0, 1) +
  labs(x = "Stimulus", y = "Probability of response", title = "Psychometric Curve") +
  theme_minimal()