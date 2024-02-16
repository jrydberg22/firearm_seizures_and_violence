#### Examining The Effect of Firearm Seizures During Directed Patrol on Firearm Violence
#### Manuscript Analysis Replication Script

#### Setup Workspace ####

# Packages 

library(tidyverse)
library(brms)
library(marginaleffects)

# Working Directory

setwd(FILEPATH)

#### Read in Data ####

d = read_csv("Firearm Seizures Replication Data.csv")

#### Setup for Post-Intervention Data ####

d = d %>%
  mutate(block_id = GEOID10,
         intv_area = hotspot)

d_post = d %>%
  filter(time > 24) %>%
  mutate(time = time - (min(time)-1))

### Setup for Reverse Time Order Checks
# rto = "Reverse Time Order"

d_pre = d %>%
  filter(time < 25) %>%
  rename(rto_total_gun = total_gun,
         rto_total_nongun = total_nongun,
         rto_gun_homicides = gun_homicides,
         rto_nongun_homicides = nongun_homicides,
         rto_gun_assaults = gun_assaults,
         rto_nongun_assaults = nongun_assaults,
         rto_gun_robberies = gun_robberies,
         rto_nongun_robberies = nongun_robberies,
         rto_weapons_offense = weapons_offense,
         rto_cfs_domestic = cfs_domestic,
         rto_cfs_shotsfired = cfs_shotsfired) %>%
  select(block_id, time, starts_with("rto_"))

# Merge to Post-Intv Data

d_post = d_post %>%
  left_join(d_pre)

#### Descriptive Statistics ####

desc_stats = function(x) {
  n = sum(!is.na(x))
  mn = mean(x, na.rm = TRUE)
  sd = sd(x, na.rm = TRUE)
  min = min(x, na.rm = TRUE)
  max = max(x, na.rm = TRUE)
  df = tibble(
    n = n,
    mean = mn,
    std_dev = sd,
    min = min,
    max = max
  )
  return(df)
}

#### Table 1
# Dependent Variables
desc_stats(d_post$total_gun)
desc_stats(d_post$cfs_shotsfired)

# Alternative Dependent Variables
desc_stats(d_post$total_nongun)
desc_stats(d_post$cfs_domestic)

# Reverse Time Order Check Dependent Variables
desc_stats(d_post$rto_total_gun)
desc_stats(d_post$rto_cfs_shotsfired)

# Intervention Area Specific Measures
htspt = d_post %>%
  distinct(intv_area, time, tx, gun_seizures, traffic_sqmi) 

desc_stats(htspt$gun_seizures)

htspt %>%
  group_by(tx) %>%
  summarize(total = sum(gun_seizures),
            mean = mean(gun_seizures))

# Checks on Firearm Seizures in Final Quarter
htspt %>%
  filter(time > 21) %>%
  group_by(intv_area) %>%
  summarize(total = sum(gun_seizures),
            mean = mean(gun_seizures))

# Control Variables
desc_stats(htspt$traffic_sqmi)

blocks = d_post %>%
  distinct(block_id, tx, t.pop, t.hu, t.hh, t.fam,
           n.femfam, n.vacant, n.renter, n.black, 
           n.hispanic, n.male1521)

desc_stats(blocks$tx)
desc_stats(blocks$t.pop)
desc_stats(blocks$n.femfam)
desc_stats(blocks$n.vacant)
desc_stats(blocks$n.black)
desc_stats(blocks$n.hispanic)
desc_stats(blocks$n.male1521)

#### Fake Data Simulation ####

library(MASS)

# Set seed for reproducibility
set.seed(123)

# Define number of neighborhoods and census blocks
n_neighborhoods = 10
n_blocks = 3000

# Define number of months and create a sequence
n_months = 24
months = 1:n_months

# Create neighborhood-level data
neighborhoods = tibble(
  neighborhood_id = 1:n_neighborhoods) %>%
  tidyr::expand(neighborhood_id, month = months) %>%
  mutate(firearms_seized = rnbinom(n = n(), mu = 1, size = 0.1))

# Create block-level data
blocks = tibble(
  block_id = 1:n_blocks,
  neighborhood_id = sample(1:n_neighborhoods, size = n_blocks, replace = TRUE),
  block_population = rnorm(n_blocks, 34, 44),
  fem_hh = rnorm(n_blocks, 4, 7)
)

# Create monthly crime data for each block
crime = blocks %>%
  tidyr::expand(block_id, month = months) %>%
  left_join(blocks, by = "block_id") %>%
  left_join(neighborhoods) %>%
  mutate(crimes = rnbinom(n = n(), mu = 1 + (0.25*firearms_seized), size = 0.01))

# Preview the data
head(crime)

rm(blocks, neighborhoods, months, n_blocks, n_months, n_neighborhoods)

### Modeling Fake Data

crime = crime %>%
  group_by(neighborhood_id) %>%
  mutate(fs_between = mean(firearms_seized)) %>%
  ungroup() %>%
  mutate(fs_within = firearms_seized - fs_between,
         monthZ = robustHD::standardize(month),
         popZ = robustHD::standardize(block_population),
         femhhZ = robustHD::standardize(fem_hh),
         obs_id = 1:nrow(crime))

# Negative Binomial

f1 = brm(crimes ~ monthZ + fs_within + fs_between + popZ + 
           (1 + monthZ | neighborhood_id) + (1 + monthZ | block_id),
         data = crime,
         family = negbinomial(),
         prior = c(prior(normal(0, 1), class = Intercept),
                   prior(normal(0, 0.5), class = b),
                   prior(exponential(0.02), class = shape),
                   prior(exponential(2), class = sd),
                   prior(lkj(4), class = cor)),
         cores = parallel::detectCores(),
         chains = 4, iter = 2000, warmup = 1000, threads = threading(2),
         control = list(adapt_delta = 0.99,
                        max_treedepth = 15),
         save_pars = save_pars(),
         backend = "cmdstanr",
         file = "fake data sim - negbin",
         file_refit = "on_change")

print(f1, digits = 3)

# Zero-Inflated RE Poisson

f2 = brm(bf(crimes ~ monthZ + fs_within + fs_between + popZ + 
              (1 + monthZ | neighborhood_id) + (1 + monthZ | block_id),
            zi ~ 1 + (1 | neighborhood_id) + (1 | block_id)),
         data = crime,
         family = zero_inflated_poisson(),
         prior = c(prior(normal(0, 1), class = Intercept),
                   prior(normal(0, 0.5), class = b),
                   prior(exponential(2), class = sd),
                   prior(lkj(4), class = cor),
                   prior(logistic(0, 1), class = Intercept, dpar = zi),
                   prior(exponential(2), class = sd, dpar = zi)),
         cores = parallel::detectCores(),
         chains = 4, iter = 2000, warmup = 1000,
         save_pars = save_pars(),
         backend = "cmdstanr",
         file = "fake data sim - zi poisson",
         file_refit = "on_change")

print(f2, digits = 3)

(waic_nb = waic(f1))
(waic_zip = waic(f2))

loo_compare(waic_nb, waic_zip)

# Compare

(loo_nb = loo(f1))
(loo_orep = loo(f2))

loo_compare(loo_nb, loo_orep)

# Zero Inflated Negbin

f3 = brm(bf(crimes ~ monthZ + fs_within + fs_between + popZ + 
              (1 + monthZ | neighborhood_id) + (1 + monthZ | block_id),
            zi ~ 1 + (1 | neighborhood_id) + (1 | block_id)),
         data = crime,
         family = zero_inflated_negbinomial(),
         prior = c(prior(normal(0, 1), class = Intercept),
                   prior(normal(0, 0.5), class = b),
                   prior(exponential(2), class = sd),
                   prior(exponential(0.02), class = shape),
                   prior(lkj(4), class = cor),
                   prior(logistic(0, 1), class = Intercept, dpar = zi),
                   prior(exponential(2), class = sd, dpar = zi)),
         cores = parallel::detectCores(),
         chains = 4, iter = 2000, warmup = 1000, threads = 2,
         save_pars = save_pars(),
         backend = "cmdstanr",
         file = "fake data sim - zi negbin",
         file_refit = "on_change")

(waic_nb = waic(f1))
(waic_zinb = waic(f3))

loo_compare(waic_nb, waic_zinb)

#### Prior Predictive Simulation ####

d_post = d_post %>%
  group_by(intv_area) %>%
  mutate(fs_between = mean(gun_seizures)) %>%
  ungroup() %>%
  mutate(fs_within = gun_seizures - fs_between,
         timeZ = robustHD::standardize(time),
         traffic_sqmiZ = robustHD::standardize(traffic_sqmi),
         t.popZ = robustHD::standardize(t.pop),
         femhhZ = robustHD::standardize(n.femfam),
         vacantZ = robustHD::standardize(n.vacant),
         blackZ = robustHD::standardize(n.black),
         hispanicZ = robustHD::standardize(n.hispanic),
         ml1521Z = robustHD::standardize(n.male1521),
         obs_id = 1:nrow(d_post))

prior = c(prior(exponential(0.02), class = shape),
          prior(normal(0, 1), class = Intercept),
          prior(normal(0, 0.1), class = b),
          prior(exponential(2), class = sd),
          prior(lkj(4), class = cor))

p1 = brm(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
           t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
           (1 + timeZ | intv_area) + (1 + timeZ | block_id),
         data = d_post,
         family = negbinomial(link = "log"),
         prior = prior,
         cores = parallel::detectCores(),
         chains = 4, iter = 2000, warmup = 1000, threads = 2,
         seed = 483313,
         backend = "cmdstanr",
         sample_prior = "only")

p1 %>%
  sjPlot::plot_model(type = "pred", terms = "fs_within") +
  theme_classic()

rm(p1)

#### Functional Form of Time ####

# Gun Violence

linear = brm(total_gun ~ timeZ +  
               (1 | intv_area) + (1 | block_id),
             data = d_post,
             family = negbinomial(link = "log"),
             prior = c(prior(exponential(0.02), class = shape),
                       prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd)),
             cores = parallel::detectCores(),
             chains = 4, iter = 2000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 299583,
             backend = "cmdstanr")

quad = brm(total_gun ~ timeZ + I(timeZ^2) + 
             (1 | intv_area) + (1 | block_id),
           data = d_post,
           family = negbinomial(link = "log"),
           prior = c(prior(exponential(0.02), class = shape),
                     prior(normal(0, 1), class = Intercept),
                     prior(normal(0, 0.1), class = b),
                     prior(exponential(2), class = sd)),
           cores = parallel::detectCores(),
           chains = 4, iter = 2000, warmup = 1000, threads = 2,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 15),
           seed = 596692,
           backend = "cmdstanr")

lin_random = brm(total_gun ~ timeZ +  
                   (1 + timeZ | intv_area) + (1 + timeZ | block_id),
                 data = d_post,
                 family = negbinomial(link = "log"),
                 prior = c(prior(exponential(0.02), class = shape),
                           prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(lkj(4), class = cor)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 181351,
                 backend = "cmdstanr")

quad_random = brm(total_gun ~ timeZ + I(timeZ^2) + 
                    (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                  data = d_post,
                  family = negbinomial(link = "log"),
                  prior = c(prior(exponential(0.02), class = shape),
                            prior(normal(0, 1), class = Intercept),
                            prior(normal(0, 0.1), class = b),
                            prior(exponential(2), class = sd),
                            prior(lkj(4), class = cor)),
                  cores = parallel::detectCores(),
                  chains = 4, iter = 2000, warmup = 1000, threads = 2,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15),
                  seed = 407013,
                  backend = "cmdstanr")

# Compare

(waic_linear = waic(ear))
(waic_quad = waic(quad))
(waic_linran = waic(lin_random))
(waic_linrwquad = waic(lin_random_wquad))
(waic_quadran = waic(quad_random))

loo_compare(waic_linear, waic_quad, waic_linran, waic_quadran)

# Shots fired

linear = brm(cfs_shotsfired ~ timeZ +  
               (1 | intv_area) + (1 | block_id),
             data = d_post,
             family = negbinomial(link = "log"),
             prior = c(prior(exponential(0.02), class = shape),
                       prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd)),
             cores = parallel::detectCores(),
             chains = 4, iter = 2000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 718990,
             backend = "cmdstanr")

quad = brm(cfs_shotsfired ~ timeZ + I(timeZ^2) + 
             (1 | intv_area) + (1 | block_id),
           data = d_post,
           family = negbinomial(link = "log"),
           prior = c(prior(exponential(0.02), class = shape),
                     prior(normal(0, 1), class = Intercept),
                     prior(normal(0, 0.1), class = b),
                     prior(exponential(2), class = sd)),
           cores = parallel::detectCores(),
           chains = 4, iter = 2000, warmup = 1000, threads = 2,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 15),
           seed = 488371,
           backend = "cmdstanr")

lin_random = brm(cfs_shotsfired ~ timeZ +  
                   (1 + timeZ | intv_area) + (1 + timeZ | block_id),
                 data = d_post,
                 family = negbinomial(link = "log"),
                 prior = c(prior(exponential(0.02), class = shape),
                           prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(lkj(4), class = cor)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 160420,
                 backend = "cmdstanr")

quad_random = brm(cfs_shotsfired ~ timeZ + I(timeZ^2) + 
                    (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                  data = d_post,
                  family = negbinomial(link = "log"),
                  prior = c(prior(exponential(0.02), class = shape),
                            prior(normal(0, 1), class = Intercept),
                            prior(normal(0, 0.1), class = b),
                            prior(exponential(2), class = sd),
                            prior(lkj(4), class = cor)),
                  cores = parallel::detectCores(),
                  chains = 4, iter = 2000, warmup = 1000, threads = 2,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15),
                  seed = 799866,
                  backend = "cmdstanr")

# Compare
# Table 3 Statistics

(waic_linear = waic(linear))
(waic_quad = waic(quad))
(waic_linran = waic(lin_random))
(waic_quadran = waic(quad_random))

loo_compare(waic_linear, waic_quad, waic_linran, waic_quadran)

#### Modeling for Inference ####

# Table 4 Models
# Gun Violence

gv_fit = brm(total_gun ~ 
               timeZ + I(timeZ^2) + fs_between + fs_within +
               tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
               (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
             data = d_post,
             family = negbinomial(link = "log"),
             prior = c(prior(exponential(0.02), class = shape),
                       prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(lkj(4), class = cor)),
             cores = parallel::detectCores(),
             chains = 4, iter = 4000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 111992,
             file = "model fits/gun violence fit",
             file_refit = "on_change",
             backend = "cmdstanr")

#gv_fit = read_rds("model fits/gun violence fit.rds")

print(gv_fit, digits = 3)

# Calls for Service - Shots Fired

sf_fit = brm(cfs_shotsfired ~ 
               timeZ + I(timeZ^2) + fs_between + fs_within +
               tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
               (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
             data = na.omit(d_post),
             family = negbinomial(link = "log"),
             prior = c(prior(exponential(0.02), class = shape),
                       prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(lkj(4), class = cor)),
             cores = parallel::detectCores(),
             chains = 4, iter = 4000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 95069,
             file = "model fits/shots fired fit",
             file_refit = "on_change",
             backend = "cmdstanr")

print(sf_fit, digits = 3)

### Alternate Dependent Variables

ngv_fit = brm(total_nongun ~ 
               timeZ + I(timeZ^2) + fs_between + fs_within +
               tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
               (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
             data = d_post,
             family = negbinomial(link = "log"),
             prior = c(prior(exponential(0.02), class = shape),
                       prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(lkj(4), class = cor)),
             cores = parallel::detectCores(),
             chains = 4, iter = 4000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 965336,
             file = "model fits/alt dv - non gun violence fit",
             file_refit = "on_change",
             backend = "cmdstanr")

#ngv_fit = read_rds("model fits/alt dv - non gun violence fit.rds")

print(ngv_fit, digits = 3)

dcfs_fit = brm(cfs_domestic ~ 
                 timeZ + I(timeZ^2) + fs_between + fs_within +
                 tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                 (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
               data = na.omit(d_post),
               family = negbinomial(link = "log"),
               prior = c(prior(exponential(0.02), class = shape),
                         prior(normal(0, 1), class = Intercept),
                         prior(normal(0, 0.1), class = b),
                         prior(exponential(2), class = sd),
                         prior(lkj(4), class = cor)),
               cores = parallel::detectCores(),
               chains = 4, iter = 4000, warmup = 1000, threads = 2,
               control = list(adapt_delta = 0.99,
                              max_treedepth = 15),
               seed = 547622,
               file = "model fits/alt dv - domestic calls fit",
               file_refit = "on_change",
               backend = "cmdstanr")

print(dcfs_fit, digits = 3)

# Reverse Time Order - Gun Violence

rtogv_fit = brm(rto_total_gun ~ 
                  timeZ + I(timeZ^2) + fs_between + fs_within +
                  tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                  (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                data = d_post,
                family = negbinomial(link = "log"),
                prior = c(prior(exponential(0.02), class = shape),
                          prior(normal(0, 1), class = Intercept),
                          prior(normal(0, 0.1), class = b),
                          prior(exponential(2), class = sd),
                          prior(lkj(4), class = cor)),
                cores = parallel::detectCores(),
                chains = 4, iter = 4000, warmup = 1000, threads = 2,
                control = list(adapt_delta = 0.99,
                               max_treedepth = 15),
                seed = 493006,
                file = "model fits/rto - gun violence fit",
                file_refit = "on_change",
                backend = "cmdstanr")

print(rtogv_fit, digits = 3)

rtocfs_fit = brm(rto_cfs_shotsfired ~ 
                   timeZ + I(timeZ^2) + fs_between + fs_within +
                   tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                   (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                 data = na.omit(d_post),
                 family = negbinomial(link = "log"),
                 prior = c(prior(exponential(0.02), class = shape),
                           prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(lkj(4), class = cor)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 4000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 965336,
                 file = "model fits/rto - shotsfired calls fit",
                 file_refit = "on_change",
                 backend = "cmdstanr")

print(rtocfs_fit, digits = 3)

#### Compute Differences in Posteriors ####

library(tidybayes)

# Extract Posteriors

gv_fsw = gv_fit %>%
  spread_draws(b_fs_within) 

ngv_fsw = ngv_fit %>%
  spread_draws(b_fs_within)

rto_gv_fsw = rtogv_fit %>%
  spread_draws(b_fs_within)

## Differences in Posteriors

gv_fsw = gv_fsw %>%
  arrange(b_fs_within)

ngv_fsw = ngv_fsw %>%
  arrange(b_fs_within)

gv_comp = data.frame(
  gv = gv_fsw$b_fs_within,
  ngv = ngv_fsw$b_fs_within
)

rto_comp = data.frame(
  gv = gv_fsw$b_fs_within,
  rto = rto_gv_fsw$b_fs_within
)

hypothesis(gv_comp, "gv = ngv", alpha = 0.05)

hypothesis(rto_comp, "gv = rto", alpha = 0.05)

plot(hypothesis(gv_fsw, "b_fs_within < -0.007", alpha = 0.05))

(gv_fsw$b_fs_within - ngv_fsw$b_fs_within) %>%
  mean_qi(.width = 0.95)

(gv_fsw$b_fs_within - rto_gv_fsw$b_fs_within) %>%
  mean_qi(.width = 0.95)

# Visualize

tibble(
  diff = c(gv_fsw$b_fs_within - ngv_fsw$b_fs_within,
            gv_fsw$b_fs_within - rto_gv_fsw$b_fs_within),
  outcome = rep(c("Alternate\nDependent\nVariable\nNon-Gun Violence", 
                  "Reverse\nTime Order\nGun Violence"), each = 12000)
) %>%
  ggplot(aes(y = outcome, x = diff)) +
  stat_halfeye(.width = c(0.9, 0.95)) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = expression(atop("Difference in Posterior Distributions", paste("Main Result ", beta[w], " - Comparison ", beta[w]))),
       y = "Robustness Check Comparison") +
  theme_classic() +
  theme(text = element_text(color = "black", family = "serif"),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16))

#### Check with Lagged Seizures ####

d_post = d_post %>%
  group_by(intv_area) %>%
  mutate(lag_seizures = lag(gun_seizures),
         fs_between_lag = mean(lag_seizures, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(fs_within_lag = lag_seizures - fs_between_lag)

d_post %>%
  select(time, block_id, gun_seizures, lag_seizures, fs_between_lag, fs_within_lag) %>%
  head(20)

## Gun Violence

gv_fit_lag = brm(total_gun ~ 
                   timeZ + I(timeZ^2) + fs_between_lag + fs_within_lag +
                   tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                   (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                 data = d_post,
                 family = negbinomial(link = "log"),
                 prior = c(prior(exponential(0.02), class = shape),
                           prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(lkj(4), class = cor)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 4000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 854332,
                 file = "model fits/gun violence fit lag",
                 file_refit = "on_change",
                 backend = "cmdstanr")

#gv_fit_lag = read_rds("model fits/gun violence fit lag.rds")

print(gv_fit_lag, digits = 3)

# Calls for Service - Shots Fired

sf_fit_lag = brm(cfs_shotsfired ~ 
                   timeZ + I(timeZ^2) + fs_between_lag + fs_within_lag +
                   tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                   (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                 data = na.omit(d_post),
                 family = negbinomial(link = "log"),
                 prior = c(prior(exponential(0.02), class = shape),
                           prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(lkj(4), class = cor)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 4000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 516009,
                 file = "model fits/shots fired fit lag",
                 file_refit = "on_change",
                 backend = "cmdstanr", silent = 0)

print(sf_fit_lag, digits = 3)

#### Cross Check with ML Poisson
# Warning: These take just as long, if not longer, than the Bayesian
# negative binomial models

library(sandwich)
library(clubSandwich)
library(lmtest)

# Gun Violence
pois_fit = glm(total_gun ~ gun_seizures + factor(time) + factor(block_id),
               data = d_post,
               family = poisson(link = "log"))
summary(pois_fit)
coeftest(pois_fit, vcov = vcovCL, cluster = ~ time + intv_area)
# -0.032 (.009) z = -3.61, p = .0003

# Gun Violence with block trends
pois_fit2 = glm(total_gun ~ gun_seizures + factor(time) + factor(block_id) +
                 (factor(block_id)*time),
               data = d_post,
               family = poisson(link = "log"))
summary(pois_fit2)
coeftest(pois_fit2, vcov = vcovCL, cluster = ~ time + intv_area)
coefci(pois_fit2, vcov = vcovCL, cluster = ~ time + block_id + intv_area,
       parm = "gun_seizures")
# # -0.031 (.008) z = -3.85, p = .0001 [-0.046, -0.015]

# Shots fired
pois_fit_cfs = glm(cfs_shotsfired ~ gun_seizures + factor(time) + factor(block_id),
               data = d_post,
               family = poisson(link = "log"))
summary(pois_fit_cfs)
coeftest(pois_fit_cfs, vcov = vcovCL, cluster = ~ time + intv_area)
# -0.001 (.003) z = -0.266, p = .791

# Shots fired with block trends
pois_fit_cfs2 = glm(cfs_shotsfired ~ gun_seizures + factor(time) + factor(block_id) +
                      (factor(block_id)*time),
                   data = d_post,
                   family = poisson(link = "log"))
summary(pois_fit_cfs2)
coeftest(pois_fit_cfs2, vcov = vcovCL, cluster = ~ time + block_id + intv_area)
coefci(pois_fit_cfs2, vcov = vcovCL, cluster = ~ time + block_id + intv_area,
       parm = "gun_seizures")
# 0.019 (.026) z = .719, p = .472

### Sensitivity of effect of reported gun violence to exclusion of final quarter

gv_fit_no13Q4 = brm(total_gun ~ 
                      timeZ + I(timeZ^2) + fs_between + fs_within +
                      tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                      (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                    data = d_post %>% filter(time < 22),
                    family = negbinomial(link = "log"),
                    prior = c(prior(exponential(0.02), class = shape),
                              prior(normal(0, 1), class = Intercept),
                              prior(normal(0, 0.1), class = b),
                              prior(exponential(2), class = sd),
                              prior(lkj(4), class = cor)),
                    cores = parallel::detectCores(),
                    chains = 4, iter = 4000, warmup = 1000, threads = 2,
                    control = list(adapt_delta = 0.99,
                                   max_treedepth = 15),
                    seed = 287578,
                    file = "model fits/gun violence fit no 2013 Q4",
                    file_refit = "on_change",
                    backend = "cmdstanr")

print(gv_fit_no13Q4, digits = 3)
