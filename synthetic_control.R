library(tidyverse)
library()
library(cmdstanr)
library(brms)
library(posterior)
library(haven)
library(Synth)
library(devtools)
if(!require(SCtools)) devtools::install_github("bcastanho/SCtools")
library(SCtools)


# function for quantiles
get_qs <- function(draws){
  
  tmp <- apply(draws, 2, function(x) quantile(x, c(0.025, .5, 0.975))) * adj_factor
  tmp <- t(tmp) 
  tmp <- data.frame(tmp)
  tmp$true <- rep(NA, nrow(tmp))
  tmp$true <- c(y_train$statefips_48, y_test$statefips_48) * adj_factor
  tmp$years <- c(1985:2000)
  
  return(tmp)
  }

read_data <- function(df)
{
  full_path <- paste("https://raw.github.com/scunning1975/mixtape/master/", 
                     df, sep = "")
  df <- read_dta(full_path)
  return(df)
}

texas <- read_data("texas.dta") %>%
  as.data.frame(.)

dataprep_out <- dataprep(
  foo = texas,
  predictors = c("poverty", "income"),
  predictors.op = "mean",
  time.predictors.prior = 1985:1993,
  special.predictors = list(
    list("bmprison", c(1988, 1990:1992), "mean"),
    list("alcohol", 1990, "mean"),
    list("aidscapita", 1990:1991, "mean"),
    list("black", 1990:1992, "mean"),
    list("perc1519", 1990, "mean")),
  dependent = "bmprison",
  unit.variable = "statefip",
  unit.names.variable = "state",
  time.variable = "year",
  treatment.identifier = 48,
  controls.identifier = c(1,2,4:6,8:13,15:42,44:47,49:51,53:56),
  time.optimize.ssr = 1985:1993,
  time.plot = 1985:2000
)

synth_out <- synth(data.prep.obj = dataprep_out)

path.plot(synth_out, dataprep_out)


gaps.plot(synth_out, dataprep_out)


placebos <- generate.placebos(dataprep_out, synth_out, Sigf.ipop = 3)

plot_placebos(placebos)

mspe.plot(placebos, discard.extreme = TRUE, mspe.limit = 1, plot.hist = TRUE)





# try BSCM is Stan.
# Prep data.
control_fips <- c(1,2,4:6,8:13,15:42,44:51,53:56)
adj_factor <- 10000
texas$bmpercap <- texas$bmprison / texas$bmpop
texas_bscm <- texas[texas$statefip %in% control_fips, c("year", "statefip", "bmprison")] %>%
  pivot_wider(names_from= statefip, values_from = bmprison,
              names_prefix = "statefips_")
X <- texas_bscm[,str_detect(colnames(texas_bscm), "statefips_48") == FALSE]
y <- texas_bscm[, "statefips_48"]
X_train <- X[X$year <= 1993, 2:51] / adj_factor
X_test <- X[X$year > 1993, 2:51] / adj_factor
y_train <- y[1:9, ] / adj_factor
y_test <- y[10:16, ] / adj_factor


bscm_data <- list(
  
  N_train = nrow(X_train),
  N_test = nrow(X_test),
  p = ncol(X_train),
  y_train = as.vector(as.matrix(y_train)),
  X_train = as.matrix(X_train),
  X_test = as.matrix(X_test)
  
)

bscm <- cmdstan_model("bscm_horseshoe.stan")

fit <- bscm$sample(
  data = bscm_data,
  seed = 321205,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 1000,
  refresh = 100,
  adapt_delta = 0.999,
  max_treedepth = 15
  
)

fit$cmdstan_diagnose()
fit$save_object(file = "fit_no_predictors.RDS")

y_test_draws <- as_draws_df(fit$draws(c("y_test", "y_fit")))
y_test_draws <- data.frame(y_test_draws[, 1:16])
#for(i in 1:7){
#  y_test_draws[, i] <- y_test_draws[, i] - y_test[i, 1]
#}

qs <- get_qs(y_test_draws)

plt <- ggplot(qs) + geom_ribbon(mapping = aes(x = years,
                                           ymin = X2.5.,
                                           ymax = X97.5.),
                             fill = "blue", 
                         alpha = 0.3) +
  geom_line(mapping = aes(x = years, 
                          y = X50.),
            colour = "blue",
            size = 2) +
  geom_line(mapping = aes(x = years, 
                          y = true),
            linetype = 2,
            size = 2) +
  geom_text(mapping = aes(x = 1996, y = 60000, label = "Texas"), size = 6) + 
  geom_text(mapping = aes(x = 1996.5, y = 40000, label = "Synthetic\nTexas"), size = 6) +
  labs(title = "Testing bayesian synthetic control",
       x = "Year",
       y = "Black male prison population") +
  geom_vline(xintercept = 1993, colour = "red", linetype = 3) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 18),
        plot.background = element_rect(fill = "white"))
ggsave(filename = "Synthetic_Texas_1.png", plot = plt,
       height = 1500, width = 2000, units = "px")
plt



# fit with new Stan model
bscm_data2 <- list(
  
  N_pre = nrow(X_train),
  N_post = nrow(X_test),
  p = ncol(X_train),
  y_pre = as.vector(as.matrix(y_train)),
  X_pre = as.matrix(X_train),
  X_post = as.matrix(X_test),
  
  scale_alpha = 1,
  scale_g = 1,
  
  slab_scale = 1,
  slab_df = 5,
  
  nu_g = 3,
  nu_l = 3
  
)

bscm2 <- cmdstan_model("bscm_horseshoe_modified.stan")

fit2 <- bscm2$sample(
  data = bscm_data2,
  seed = 321205,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 100,
  adapt_delta = 0.9999,
  max_treedepth = 15
  
)

fit2$cmdstan_diagnose()
fit2$save_object(file = "fit2_no_predictors.RDS")



y_test_draws2 <- as_draws_df(fit2$draws(c("y_fit", "y_post")))
y_test_draws2 <- y_difference_draws <- data.frame(y_test_draws[, 1:16])


y_all <- rbind(y_train, y_test)
for(i in 1:16){
  y_difference_draws[, i] <- y_all[i, 1] - y_difference_draws[, i]
}

qs2 <- get_qs(y_test_draws2)

plt2 <- ggplot(qs2) + geom_ribbon(mapping = aes(x = years,
                                              ymin = X2.5.,
                                              ymax = X97.5.),
                                fill = "blue", 
                                alpha = 0.3) +
  geom_line(mapping = aes(x = years, 
                          y = X50.),
            colour = "blue",
            size = 2) +
  geom_line(mapping = aes(x = years, 
                          y = true),
            linetype = 2,
            size = 2) +
  geom_text(mapping = aes(x = 1996, y = 60000, label = "Texas"), size = 6) + 
  geom_text(mapping = aes(x = 1996.5, y = 40000, label = "Synthetic\nTexas"), size = 6) +
  labs(title = "Testing bayesian synthetic control",
       subtitle = "Stan script",
       x = "Year",
       y = "Black male prison population") +
  geom_vline(xintercept = 1993, colour = "red", linetype = 3) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        plot.background = element_rect(fill = "white"))
ggsave(filename = "Synthetic_Texas_2.png", plot = plt2,
       height = 1500, width = 2000, units = "px")
plt2



qs3 <- get_qs(y_difference_draws)

plt3 <- ggplot(qs3) + geom_ribbon(mapping = aes(x = years,
                                                ymin = X2.5.,
                                                ymax = X97.5.),
                                  fill = "blue", 
                                  alpha = 0.3) +
  geom_line(mapping = aes(x = years, 
                          y = X50.),
            colour = "blue",
            size = 2) +
  labs(title = "Testing bayesian synthetic control",
       subtitle = "Difference from predicted trend",
       x = "Year",
       y = "Black male prison population") +
  geom_vline(xintercept = 1993, colour = "red", linetype = 3) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        plot.background = element_rect(fill = "white"))
ggsave(filename = "Synthetic_Texas_3.png", plot = plt3,
       height = 1500, width = 2000, units = "px")
plt3






# Test with brms
df_brms <- data.frame(y_train, X_train)
fit3 <-  brm(df_brms,
             formula = statefips_48 ~ 1 + ., 
             family = gaussian(),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(horseshoe(df = 4, 
                                       scale_global = 1,
                                       scale_slab = 1,
                                       df_global = 3,
                                       df_slab = 4), class = b),
                       prior(exponential(2), class = sigma)),
             chains = 4,
             cores = 4, 
             iter = 2000,
             warmup = 1000,
             backend = "cmdstanr",
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15))

y_post_brms <- posterior_epred(fit3, newdata = rbind(X_train, X_test))
pp_check()


qs_brms <- get_qs(y_post_brms)


plt_brms <- ggplot(qs_brms) + geom_ribbon(mapping = aes(x = years,
                                              ymin = X2.5.,
                                              ymax = X97.5.),
                                fill = "blue", 
                                alpha = 0.3) +
  geom_line(mapping = aes(x = years, 
                          y = X50.),
            colour = "blue",
            size = 2) +
  geom_line(mapping = aes(x = years, 
                          y = true),
            linetype = 2,
            size = 2) +
  geom_text(mapping = aes(x = 1996, y = 60000, label = "Texas"), size = 6) + 
  geom_text(mapping = aes(x = 1996.5, y = 40000, label = "Synthetic\nTexas"), size = 6) +
  labs(title = "Testing bayesian synthetic control",
       subtitle = "Using brms",
       x = "Year",
       y = "Black male prison population") +
  geom_vline(xintercept = 1993, colour = "red", linetype = 3) +
  theme_minimal() + 
  theme(plot.title = element_text(size = 18),
        plot.subtitle = element_text(size = 14),
        plot.background = element_rect(fill = "white"))

plt_brms
ggsave(filename = "Synthetic_Texas_brms.png", plot = plt_brms,
       height = 1500, width = 2000, units = "px")




