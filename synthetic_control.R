library(tidyverse)
library(cmdstanr)
library(brms)
library(posterior)
library(haven)
library(Synth)
library(devtools)
if(!require(SCtools)) devtools::install_github("bcastanho/SCtools")
library(SCtools)


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
  iter_warmup = 2000,
  iter_sampling = 1000,
  refresh = 100,
  adapt_delta = 0.999,
  max_treedepth = 15
  
)

fit$cmdstan_diagnose()
fit$save_object(file = "fit_no_predictors.RDS")

y_test_draws <- as_draws_df(fit$draws("y_test"))
y_test_draws <- data.frame(y_test_draws[, 1:7])
for(i in 1:7){
  y_test_draws[, i] <- y_test_draws[, i] - y_test[i, 1]
}

qs <- apply(y_test_draws, 2, function(x) quantile(x, c(0.025, .5, 0.975))) * adj_factor
qs <- t(qs) 
qs <- data.frame(qs)
qs$true <- y_test$statefips_48 * adj_factor
qs$years <- 1994:2000

ggplot(qs) + geom_ribbon(mapping = aes(x = years,
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
  labs(title = "Testing bayesian synthetic control",
       x = "Year",
       y = "Difference") +
  theme_minimal()
