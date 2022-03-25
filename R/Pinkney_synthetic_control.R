library(tidyverse)
library(cmdstanr)
library(data.table)
library(posterior)

load("/Users/stephenwild/Downloads/smoking.rda")

# stan_file <- "synth_penalized.stan"
#stan_file <- "synth_horseshoe_b_tau_x.stan"
stan_file <- "/Users/stephenwild/Desktop/Stats stuff/Bayesian synthetic control/bscm_improved_extended_pinkney.stan"
mod_synth <- cmdstan_model(stan_file)

smoking <- data.table(smoking)
smoking_y <- dcast(smoking, year ~ state, value.var = "cigsale")
summaries <- smoking %>%
  group_by(state) %>%
  summarize(retprice = mean(retprice, na.rm = TRUE),
            lnincome = mean(lnincome, na.rm = TRUE),
            age15to24 = mean(age15to24, na.rm = TRUE),
            beer = mean(beer, na.rm = TRUE))

sm1975 <- data.frame(smoking[year == 1975, "cigsale"],
                     smoking[year == 1975, "state"])
sm1980 <- data.frame(smoking[year == 1980, "cigsale"],
                     smoking[year == 1980, "state"])
sm1988 <- data.frame(smoking[year == 1988, "cigsale"],
                     smoking[year == 1988, "state"])
sm1975 <- rename(sm1975, cigsale1975 = cigsale)
sm1980 <- rename(sm1980, cigsale1980 = cigsale)
sm1988 <- rename(sm1988, cigsale1988 = cigsale)

summaries <- inner_join(summaries, sm1975, by = "state")
summaries <- inner_join(summaries, sm1980, by = "state")
summaries <- inner_join(summaries, sm1988, by = "state")

predictors_target <- t(as.matrix(summaries[summaries$state == "California", 2:8]))
predictors_control <- t(as.matrix(summaries[summaries$state != "California", 2:8]))
X_pred <- cbind(predictors_target, predictors_control)

target <- "California"
target_index <- which(names(smoking_y) == target)
other_index <- which(names(smoking_y) %in% names(smoking_y)[c(-1, -4)])

which(smoking_y$year == 1988)


X_pred[3, ] <- log(X_pred[3, ] / (1 - X_pred[3, ]))

stan_data <- list(
  T = nrow(smoking_y),
  J = ncol(smoking_y) - 1,
  L = 8,
  P = nrow(X_pred),
  X = X_pred,
  Y = transpose(smoking_y[ , c(..target_index, ..other_index)]),
  trt_times = nrow(smoking_y) - which(smoking_y$year == 1988)
)

fit_synth <- mod_synth$sample(data = stan_data,
                              seed = 123123,
                              iter_warmup = 500,
                              iter_sampling = 500,
                              init = 0.1,
                              chains = 4,
                              parallel_chains = 4,
                              max_treedepth = 13,
                              adapt_delta = 0.8
)


synth_out <- data.frame(as_draws_df(fit_synth$draws(c("synth_out")))[1:(1209-3)])
predictors_out <- data.frame(as_draws_df(fit_synth$draws(c("chi")))[1:7])

predictors_summary <- data.frame(names = names(summaries[,2:8]),
                                 mean = apply(predictors_out, 2, mean),
                                 lower = apply(predictors_out, 2, function(x) quantile(x, 0.025)),
                                 upper = apply(predictors_out, 2, function(x) quantile(x, 0.975)))

params <- ggplot(predictors_summary) + 
  geom_pointrange(mapping = aes(x = reorder(names, mean),
                                y = mean,
                                ymin = lower, 
                                ymax = upper)) +
  geom_vline(xintercept = 0,
             color = "red",
             linetype = 2) +
  labs(title = "Parameter estimates",
       xlab = NULL,
       ylab = "Coefficient") +
  coord_flip() +
  theme_minimal()


indx <- seq(1, 1209, by = 39)
synth_cal <- synth_out[,indx]

synth_summary <- data.frame(year = 1970:2000,
                           mean = apply(synth_cal, 2, mean),
                           lower = apply(synth_cal, 2, function(x) quantile(x, 0.05)),
                           upper = apply(synth_cal, 2, function(x) quantile(x, 0.95)),
                           l25 = apply(synth_cal, 2, function(x) quantile(x, 0.25)),
                           u75 = apply(synth_cal, 2, function(x) quantile(x, 0.75)))


synth_plot <- ggplot() + 
  geom_ribbon(data = synth_summary,
              mapping = aes(x = year,
                            ymin = lower, 
                            ymax = upper),
              fill = "blue",
              alpha = 0.5) +
  geom_pointrange(data = synth_summary,
                  mapping = aes(x = year,
                                y = mean,
                                ymin = l25, 
                                ymax = u75),
                  size = 2,
                  fatten = 0.5) +
  geom_pointrange(data = synth_summary,
                  mapping = aes(x = year,
                                y = mean,
                                ymin = lower, 
                                ymax = upper)) +
  geom_line(data = smoking_y[, c("year", "California")],
            mapping = aes(x = year,
                          y = California),
            colour = "black") +
  labs(title = "Synthetic control on California tobacco data",
       subtitle = "90% interval",
       x = "Year",
       y = "Num. packs") +
  theme_minimal()

synth_plot



read_data <- function(df)
{
  full_path <- paste("https://raw.github.com/scunning1975/mixtape/master/", 
                     df, sep = "")
  df <- read_dta(full_path)
  return(df)
}

texas <- read_data("texas.dta") %>%
  as.data.frame(.)

control_fips <- c(1,2,4:6,8:13,15:42,44:47, 49:51,53:56)
texas_fip <- 48
Y_donor <- texas[texas$statefip %in% control_fips, c("statefip", "bmprison", "year")] %>%
  pivot_wider(names_from= statefip, values_from = bmprison,
              names_prefix = "statefips_")
X_donor <- texas %>%
  group_by(statefip) %>%
  summarize(black = mean(black, na.rm = TRUE),
            alcohol = mean(alcohol, na.rm = TRUE),
            perc1519 = mean(perc1519, na.rm = TRUE),
            aidscapita = mean(aidscapita, na.rm = TRUE))
Y_target <- texas[texas$statefip == texas_fip, c("statefip", "bmprison", "year")] %>%
  pivot_wider(names_from= statefip, values_from = bmprison,
              names_prefix = "statefips_")

X_target <- X_donor[X_donor$statefip == texas_fip, 2:5]
X_donor <- X_donor[X_donor$statefip != texas_fip, 2:5]
X_all <- rbind(X_target, X_donor)


Y_all <- cbind(Y_target[,2], Y_donor[, 2:51])
X_all <- t(as.matrix(X_all))


bscm_data3 <- list(
  
  T = nrow(Y_all),
  J = ncol(Y_all),
  L = 10,
  P = nrow(X_all),
  X = X_all,
  Y = as.matrix(transpose(Y_all / 10000)),
  trt_times = 7
  
)

fit_synth_texas <- mod_synth$sample(data = bscm_data3,
                              seed = 123123,
                              iter_warmup = 500,
                              iter_sampling = 500,
                              init = 0.1,
                              chains = 4,
                              parallel_chains = 4,
                              max_treedepth = 13,
                              adapt_delta = 0.8
)



install.packages("gapminder")
library(gapminder)
mod2 <- brm(lifeExp ~ 1 + year + (1 + year | country), data = dfg, cores = 4, chains = 4, backend = "cmdstanr", control = list(adapt_delta = 0.95))




