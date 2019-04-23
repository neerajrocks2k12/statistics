library("ggplot2")
library(dplyr)

options(scipen = 999)

head(faithful)
help("faithful")

p <- ggplot(faithful, aes(x = waiting)) + geom_density()
p

p + geom_vline(xintercept = 53, col = "red", size = 1) +
    geom_vline(xintercept = 80, col = "blue", size = 1)
install.packages("mixtools")
library(mixtools)

#' @param x Input data
#' @param mu mean of the component
#' @param sigma is std dev of the component
#' @param lam is lambda or the weight of the component
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

set.seed(1)
wait <- faithful$waiting
mix.Modal <- normalmixEM(wait, k = 2)

data.frame(x = mix.Modal$x) %>%
  ggplot() + 
  geom_histogram(aes(x, ..density..), binwidth = 1, color = "blue", fill = "white") +
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix.Modal$mu[1], mix.Modal$sigma[1], mix.Modal$lambda[1]),
                colour = "red", lwd = 1.5) + 
  stat_function(geom = "line", fun = plot_mix_comps,
                args = list(mix.Modal$mu[2], mix.Modal$sigma[2], mix.Modal$lambda[2]),
                colour = "blue", lwd = 1.5) + 
  ylab("Density")

# posterior probability of belonging to either component
post.df <- as.data.frame(cbind(x = mix.Modal$x, mix.Modal$posterior))
head(post.df, 10)

post.df %>%
  filter(x > 66, x < 68)

post.df %>%
  mutate(label = ifelse(comp.1 > 0.3, 1, 2)) %>%
  ggplot(aes(x = factor(label))) +
  geom_bar() + 
  xlab("component") + 
  ylab("# of data points")

# we can use any dist not necessarily Gaussian with EM. Binomials, Multinomials, 
# student-t etc.
# EM is a soft label algo unlike k-means where the data point gets "assigned"
# to a cluster. In EM we get probability of belonging to either "component"
