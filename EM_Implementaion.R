library(ggplot2)
library(dplyr)
library(reshape2)

options(scipen = 999)
set.seed(1)

comp1.vals <- data_frame(comp = "A",
                         vals = rnorm(50, mean = 1, sd = 0.5))

comp2.vals <- data_frame(comp = "B",
                         vals = rnorm(50, mean = 1.5, sd = 0.5))

vals.df <- bind_rows(comp1.vals, comp2.vals)

vals.df %>%
  ggplot(aes(x = vals, y = "A", color = factor(comp))) +
  geom_point(size = 2, alpha = 0.4) +
  scale_color_discrete(name = "Source of Data") + 
  xlab("Values") + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "top")

vals.df %>%
  group_by(comp) %>%
  summarize(mean_vals = mean(vals),
            sd_vals = sd(vals))
  
vals.df %>%
  ggplot(aes(x = vals, y = 0)) +
  geom_point(size = 2, alpha = 0.4) + 
  xlab("values") + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

#---------------
wait <- faithful$waiting

wait.kmeans <- kmeans(wait, 2)
wait.kmeans.cluster <- wait.kmeans$cluster

wait.df <- data_frame(x = wait, cluster = wait.kmeans.cluster)
wait.df

wait.df %>%
  mutate(num = row_number()) %>%
  ggplot(aes(y = num, x = x, color = factor(cluster))) + 
  geom_point() + 
  xlab("Values") +
  ylab("Data point number") +
  scale_color_discrete(name = "cluster") +
  ggtitle("kmeans plot")

wait.summary.df <- wait.df %>%
  group_by(cluster) %>%
  summarize(mu = mean(x), variance = var(x),
            std = sd(x), size = n())
wait.summary.df %>%
  select(cluster, mu, variance, std)

wait.summary.df <- wait.summary.df %>%
  mutate(alpha = size / sum(size))

wait.summary.df %>%
  select(cluster, size, alpha)

comp1.prod <- dnorm(66, wait.summary.df$mu[1], wait.summary.df$std[1]) * 
  wait.summary.df$alpha[1]

comp2.prod <- dnorm(66, wait.summary.df$mu[2], wait.summary.df$std[2]) * 
  wait.summary.df$alpha[2]

normalizer = comp1.prod + comp2.prod
comp1.prod / normalizer

comp1.prod.all <- dnorm(x = wait, wait.summary.df$mu[1], wait.summary.df$std[1]) * 
  wait.summary.df$alpha[1]

comp2.prod.all <- dnorm(x = wait, wait.summary.df$mu[2], wait.summary.df$std[2]) * 
  wait.summary.df$alpha[2]

normalizer.all <- comp1.prod.all + comp2.prod.all

comp1.post <- comp1.prod.all / normalizer.all
comp2.post <- comp2.prod.all / normalizer.all

comp1.n <- sum(comp1.post)
comp2.n <- sum(comp2.post)

comp1.mu <- (1 / comp1.n) * sum(comp1.post * wait)
comp2.mu <- (1 / comp2.n) * sum(comp2.post * wait)

comp1.var <- (1 / comp1.n) * sum(comp1.post * (wait - comp1.mu)^2)
comp2.var <- (1 / comp2.n) * sum(comp2.post * (wait - comp2.mu)^2)

comp1.alpha <- comp1.n / length(wait)
comp2.alpha <- comp2.n / length(wait)

comp.params.df <- data.frame(comp = c("comp1", "comp2"),
                             comp.mu = c(comp1.mu, comp2.mu),
                             comp.var = c(comp1.var, comp2.var),
                             comp.alpha = c(comp1.alpha, comp2.alpha),
                             comp.cal = c("self", "self"))
comp.params.df

#-------------------------------

#' Expectation step of the EM Algorithm
#' Calculates the posterior prob. (soft labels) that each component
#' has to each data point
#' 
#' @param sd.vector Vector containing standard deviation of each component
#' @param mu.vector vector containing mean of each component
#' @param alpha.vector containing the mixing weight of each component
#' @return Named list containing loglik and posterior.df
e_step <- function(x, mu.vector, sd.vector, alpha.vector) {
  comp1.prod <- dnorm(x, mu.vector[1], sd.vector[1]) * alpha.vector[1]
  comp2.prod <- dnorm(x, mu.vector[2], sd.vector[2]) * alpha.vector[2]
  sum.of.comps <- comp1.prod + comp2.prod
  comp1.post <- comp1.prod / sum.of.comps
  comp2.post <- comp2.prod / sum.of.comps
  
  sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
  sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
  
  list("loglik" = sum.of.comps.ln,
       "posterior.df" = cbind(comp1.post, comp2.post))
}

#' Maximization Step of the EM Algorithm
#'
#' Update the Component Parameters
#'
#' @param x Input data.
#' @param posterior.df Posterior probability data.frame.
#' @return Named list containing the mean (mu), variance (var), and mixing
#'   weights (alpha) for each component.
m_step <- function(x, posterior.df) {
  comp1.n <- sum(posterior.df[, 1])
  comp2.n <- sum(posterior.df[, 2])
  
  comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
  comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)
  
  comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
  comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n
  
  comp1.alpha <- comp1.n / length(x)
  comp2.alpha <- comp2.n / length(x)
  
  list("mu" = c(comp1.mu, comp2.mu),
       "var" = c(comp1.var, comp2.var),
       "alpha" = c(comp1.alpha, comp2.alpha))
}

for (i in 1:50) {
  if (i == 1) {
    # Initialization
    e.step <- e_step(wait, wait.summary.df[["mu"]], wait.summary.df[["std"]],
                     wait.summary.df[["alpha"]])
    m.step <- m_step(wait, e.step[["posterior.df"]])
    cur.loglik <- e.step[["loglik"]]
    loglik.vector <- e.step[["loglik"]]
  } else {
    # Repeat E and M steps till convergence
    e.step <- e_step(wait, m.step[["mu"]], sqrt(m.step[["var"]]), 
                     m.step[["alpha"]])
    m.step <- m_step(wait, e.step[["posterior.df"]])
    loglik.vector <- c(loglik.vector, e.step[["loglik"]])
    
    loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
    if(loglik.diff < 1e-6) {
      break
    } else {
      cur.loglik <- e.step[["loglik"]]
    }
  }
}
length(loglik.vector)
loglik.vector[1:100]

m.step

#' Plot a Mixture Component
#' 
#' @param x Input ata.
#' @param mu Mean of component.
#' @param sigma Standard of component.
#' @param lam Mixture weight of component.
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
data.frame(x = wait) %>%
  ggplot() +
  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
                 fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
                            lam = m.step$alpha[1]),
                colour = "red", lwd = 1.5) +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
                            lam = m.step$alpha[2]),
                colour = "blue", lwd = 1.5) +
  ylab("Density") +
  xlab("Values") +
  ggtitle("Final GMM Fit")
