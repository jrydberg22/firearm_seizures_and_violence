#### Using Bayesian Mixed Effect Generalized Linear Models to Evaluate Criminological Interventions: 
#### An Application to Firearm Seizures during Directed Patrol
#### Final Data Preparation and Modeling

#### Setup Workspace ####

# Packages 

library(tidyverse)
library(brms)
library(marginaleffects)
library(tidybayes)

# Working Directory

setwd("YOUR FILEPATH HERE")

#### Read in Data ####

load("Flint Firearm Seizures Analysis Data.rds")

#### Setup for Post-Intervention ####

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

# Firearm Seizures
htspt = d_post %>%
  distinct(intv_area, time, tx, gun_seizures, traffic_sqmi) 

desc_stats(htspt$gun_seizures)

htspt %>%
  group_by(tx) %>%
  summarize(total = sum(gun_seizures),
            mean = mean(gun_seizures))

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

#### Determinining Likelihood ####

## Compute gun seizure variables, standardize time and controls

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

## Gaussian

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd),
          prior(exponential(1), class = sigma))

lik1_gauss = brm(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                   t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                   (1 | intv_area) + (1 | block_id),
                 data = d_post,
                 family = gaussian(link = "identity"),
                 prior = prior,
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 seed = 739033,
                 file = "gun violence - gauss lik",
                 file_refit = "on_change",
                 backend = "rstan")

lik1_gauss = add_criterion(lik1_gauss, "waic")

## Poisson

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd))

lik2_pois = brm(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                  t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                  (1 | intv_area) + (1 | block_id),
                data = d_post,
                family = poisson(link = "log"),
                prior = prior,
                cores = parallel::detectCores(),
                chains = 4, iter = 2000, warmup = 1000, threads = 2,
                seed = 965336,
                file = "gun violence - pois lik",
                file_refit = "on_change",
                backend = "rstan")

lik2_pois = add_criterion(lik2_pois, "waic")

## Negative Binomial

prior = c(prior(exponential(1), class = shape),
          prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd))

lik3_nb = brm(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                (1 | intv_area) + (1 | block_id),
              data = d_post,
              family = negbinomial(link = "log"),
              prior = prior,
              cores = parallel::detectCores(),
              chains = 4, iter = 2000, warmup = 1000, threads = 2,
              seed = 318176,
              file = "gun violence - nb lik",
              file_refit = "on_change",
              backend = "rstan")

lik3_nb = add_criterion(lik3_nb, "waic")

## ZI Poisson

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd),
          prior(logistic(0, 1), class = Intercept, dpar = zi),
          prior(exponential(1), class = sd, dpar = zi))

lik4_zip = brm(bf(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                    t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                    (1 | intv_area) + (1 | block_id),
                  zi ~ 1 + (1 | intv_area) + (1 | block_id)),
                 data = d_post,
                 family = zero_inflated_poisson(),
                 prior = prior,
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 seed = 182372,
                 file = "gun violence - zip lik",
                 file_refit = "on_change",
                 backend = "rstan")

lik4_zip = add_criterion(lik4_zip, "waic")

## ZI Negative Binomial

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd),
          prior(exponential(1), class = shape),
          prior(logistic(0, 1), class = Intercept, dpar = zi),
          prior(exponential(1), class = sd, dpar = zi))

lik5_zinb = brm(bf(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                     t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                     (1 | intv_area) + (1 | block_id),
                   zi ~ 1 + (1 | intv_area) + (1 | block_id)),
               data = d_post,
               family = zero_inflated_negbinomial(),
               prior = prior,
               cores = parallel::detectCores(),
               chains = 4, iter = 2000, warmup = 1000, threads = 2,
               seed = 855090,
               file = "gun violence - zinb lik",
               file_refit = "on_change",
               backend = "rstan")

lik5_zinb = add_criterion(lik5_zinb, "waic")

# Print WAIC Statistics

waic(lik1_gauss)
waic(lik2_pois)
waic(lik3_nb)
waic(lik4_zip)
waic(lik5_zinb)

# Initial comparison including Gaussian

loo_compare(lik1_gauss,
            lik2_pois,
            lik3_nb,
            lik4_zip,
            lik5_zinb, criterion = "waic")

# Subsequent comparison excluding Gaussian (right column, Table 2)

loo_compare(lik2_pois,
            lik3_nb,
            lik4_zip,
            lik5_zinb, criterion = "waic")

## Gaussian

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd),
          prior(exponential(1), class = sigma))

lik1_gauss = brm(cfs_shotsfired ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                   t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                   (1 | intv_area) + (1 | block_id),
                 data = d_post,
                 family = gaussian(link = "identity"),
                 prior = prior,
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 seed = 205531,
                 file = "model fits/cfs - gauss lik",
                 file_refit = "on_change",
                 backend = "cmdstanr")

lik1_gauss = add_criterion(lik1_gauss, "waic")
lik1_gauss = read_rds("model fits/gun violence - gauss lik.rds")

## Poisson

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd))

lik2_pois = brm(cfs_shotsfired ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                  t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                  (1 | intv_area) + (1 | block_id),
                data = d_post,
                family = poisson(link = "log"),
                prior = prior,
                cores = parallel::detectCores(),
                chains = 4, iter = 2000, warmup = 1000, threads = 2,
                seed = 426568,
                file = "model fits/cfs - pois lik",
                file_refit = "on_change",
                backend = "cmdstanr")

lik2_pois = add_criterion(lik2_pois, "waic")

## Negative Binomial

prior = c(prior(exponential(1), class = shape),
          prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd))

lik3_nb = brm(cfs_shotsfired ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                (1 | intv_area) + (1 | block_id),
              data = d_post,
              family = negbinomial(link = "log"),
              prior = prior,
              cores = parallel::detectCores(),
              chains = 4, iter = 2000, warmup = 1000, threads = 2,
              seed = 586148,
              file = "model fits/cfs - nb lik",
              file_refit = "on_change",
              backend = "cmdstanr")

lik3_nb = add_criterion(lik3_nb, "waic")

## ZI Poisson

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd),
          prior(logistic(0, 1), class = Intercept, dpar = zi),
          prior(exponential(1), class = sd, dpar = zi))

lik4_zip = brm(bf(cfs_shotsfired ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                    t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                    (1 | intv_area) + (1 | block_id),
                  zi ~ 1 + (1 | intv_area) + (1 | block_id)),
               data = d_post,
               family = zero_inflated_poisson(),
               prior = prior,
               cores = parallel::detectCores(),
               chains = 4, iter = 2000, warmup = 1000, threads = 2,
               seed = 728361,
               file = "model fits/cfs - zip lik",
               file_refit = "on_change",
               backend = "cmdstanr")

lik4_zip = add_criterion(lik4_zip, "waic")

## ZI Negative Binomial

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 0.5), class = b),
          prior(exponential(1), class = sd),
          prior(exponential(1), class = shape),
          prior(logistic(0, 1), class = Intercept, dpar = zi),
          prior(exponential(1), class = sd, dpar = zi))

lik5_zinb = brm(bf(cfs_shotsfired ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
                     t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
                     (1 | intv_area) + (1 | block_id),
                   zi ~ 1 + (1 | intv_area) + (1 | block_id)),
                data = d_post,
                family = zero_inflated_negbinomial(),
                prior = prior,
                cores = parallel::detectCores(),
                chains = 4, iter = 2000, warmup = 1000, threads = 2,
                seed = 840339,
                file = "model fits/cfs - zinb lik",
                file_refit = "on_change",
                backend = "cmdstanr")

lik5_zinb = add_criterion(lik5_zinb, "waic")

waic(lik1_gauss)
waic(lik2_pois)
waic(lik3_nb)
waic(lik4_zip)
waic(lik5_zinb)

loo_compare(lik1_gauss,
            lik2_pois,
            lik3_nb,
            lik4_zip,
            lik5_zinb, criterion = "waic")

loo_compare(lik2_pois,
            lik3_nb,
            lik4_zip,
            lik5_zinb, criterion = "waic")

### Posterior Predictive Model Comparison - What is happening with Gaussian likelihood?

library(bayesplot)

# Identify outcome vector, using gun violence for this check

gv = d_post$total_gun

# 100 Posterior draws from each model

gv_rep_gaus = posterior_predict(lik1_gauss, draws = 100)
gv_rep_pois = posterior_predict(lik2_pois, draws = 100)
gv_rep_nb = posterior_predict(lik3_nb, draws = 100)
gv_rep_zip = posterior_predict(lik4_zip, draws = 100)
gv_rep_zinb = posterior_predict(lik5_zinb, draws = 100)

# Function for identifying proportions of zeros

prop_zero = function(x) mean(x == 0)

prop_zero(gv) # 98%

## Apply stat function to posterior draws

# Proportion of zeroes

gvgaus_zero = ppc_stat(gv, gv_rep_gaus, stat = "prop_zero") 
gvpois_zero = ppc_stat(gv, gv_rep_pois, stat = "prop_zero")
gvnb_zero = ppc_stat(gv, gv_rep_nb, stat = "prop_zero")
gvzip_zero = ppc_stat(gv, gv_rep_zip, stat = "prop_zero")
gvzinb_zero = ppc_stat(gv, gv_rep_zinb, stat = "prop_zero")

# Maximum value

gvgaus_max = ppc_stat(gv, gv_rep_gaus, stat = "max")
gvpois_max = ppc_stat(gv, gv_rep_pois, stat = "max")
gvnb_max = ppc_stat(gv, gv_rep_nb, stat = "max")
gvzip_max = ppc_stat(gv, gv_rep_zip, stat = "max")
gvzinb_max = ppc_stat(gv, gv_rep_zinb, stat = "max")

# Visualize a subset (just 3 for convenience, could do all 5 but that's a busy graph)

gv_ppc = tibble(
  stat = c(rep(c("Proportion of Zeros"), times = 4000*3),
           rep(c("Maximum Value"), times = 4000*3)) %>%
    factor(levels = c("Proportion of Zeros", "Maximum Value")),
  likelihood = c(rep(c("Gaussian", "Poisson", "Zero-Inflated\nNegative Binomial"),
                     each = 4000),
                 rep(c("Gaussian", "Poisson", "Zero-Inflated\nNegative Binomial"),
                     each = 4000)),
  value = c(gvgaus_zero$data$value,
            gvpois_zero$data$value,
            gvzinb_zero$data$value,
            gvgaus_max$data$value,
            gvpois_max$data$value,
            gvzinb_max$data$value),
  actual = c(rep(prop_zero(gv), times = 4000*3),
             rep(max(gv), times = 4000*3))
)

gv_ppc %>%
  ggplot(aes(x = value)) +
  geom_histogram(color = "white", fill = "deepskyblue") +
  geom_vline(aes(xintercept = actual), linetype = 1, 
             color = "black", size = 1.25, alpha = 0.5) +
  facet_wrap(vars(stat, likelihood), nrow = 2,
             scales = "free") +
  labs(x = "Posterior Summary Statistic",
       y = "Number of Posterior Draws") +
  scale_x_continuous(n.breaks = 4) +
  theme_classic() +
  theme(text = element_text(family = "serif", size = 12),
        strip.text = element_text(face = "bold"))

ggsave("Posterior Predictive.png", width = 6.5, height = 5, units = "in")

#### Prior Predictive Simulation ####

### Initial set of priors

prior = c(prior(normal(0, 5), class = Intercept),
          prior(normal(0, 1), class = b),
          prior(exponential(1), class = sd),
          prior(exponential(1), class = shape),
          prior(lkj(1), class = cor),
          prior(logistic(0, 1), class = Intercept, dpar = zi),
          prior(exponential(1), class = sd, dpar = zi))

p1 = brm(bf(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
              t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
              (1 + timeZ | intv_area) + (1 + timeZ | block_id),
            zi ~ 1 + (1 | intv_area) + (1 | block_id)),
         data = d_post,
         family = zero_inflated_negbinomial(),
         prior = prior,
         cores = parallel::detectCores(),
         chains = 4, iter = 2000, warmup = 1000, threads = 2,
         seed = 965336,
         backend = "rstan",
         sample_prior = "only") # Turning this option on ensures *only* the prior is used

# Export predictions

p1d = p1 %>%
  sjPlot::get_model_data(type = "pred", terms = "timeZ")

p1d %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_lineribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = "deepskyblue") +
  labs(x = "Time (Standardized)",
       y = "Predicted Gun Violence Incidents") +
  theme_classic() +
  theme(text = element_text(color = "black", family = "serif"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12),
        aspect.ratio = 1)

# Left panel of Figure 2

ggsave("Prior Predictive - Initial.png", 
       width = 3.25, height = 3.25, units = "in")

### Tuned Prior - made adjustments parameter class by parameter class

prior = c(prior(normal(0, 1), class = Intercept),
          prior(normal(0, 0.1), class = b),
          prior(exponential(2), class = sd),
          prior(exponential(0.02), class = shape),
          prior(lkj(4), class = cor),
          prior(logistic(0, 1), class = Intercept, dpar = zi),
          prior(exponential(2), class = sd, dpar = zi))

p2 = brm(bf(total_gun ~ timeZ + fs_within + fs_between + tx + traffic_sqmiZ + 
              t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z + 
              (1 + timeZ | intv_area) + (1 + timeZ | block_id),
            zi ~ 1 + (1 | intv_area) + (1 | block_id)),
         data = d_post,
         family = zero_inflated_negbinomial(),
         prior = prior,
         cores = parallel::detectCores(),
         chains = 4, iter = 2000, warmup = 1000, threads = 2,
         seed = 693325,
         backend = "cmdstanr",
         sample_prior = "only")

p2d = p2 %>%
  sjPlot::get_model_data(type = "pred", terms = "timeZ")

p2d %>%
  ggplot(aes(x = x, y = predicted)) +
  geom_lineribbon(aes(ymin = conf.low, ymax = conf.high),
                  fill = "deepskyblue") +
  labs(x = "Time (Standardized)",
       y = "Predicted Gun Violence Incidents") +
  theme_classic() +
  theme(text = element_text(color = "black", family = "serif"),
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 12),
        aspect.ratio = 1)

# Right panel of Figure 2

ggsave("Prior Predictive - Tuned.png", 
       width = 3.25, height = 3.25, units = "in")

#### Functional Form of Time ####

## Gun Violence

linear = brm(bf(total_gun ~ timeZ +  
               (1 | intv_area) + (1 | block_id),
               zi ~ 1 + (1 | intv_area) + (1 | block_id)),
             data = d_post,
             family = zero_inflated_negbinomial(),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(exponential(0.02), class = shape),
                       prior(logistic(0, 1), class = Intercept, dpar = zi),
                       prior(exponential(2), class = sd, dpar = zi)),
             cores = parallel::detectCores(),
             chains = 4, iter = 2000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 299583,
             file = "gvff - zinb lin",
             file_refit = "on_change",
             backend = "rstan")

linear = add_criterion(linear, "waic")

quad = brm(bf(total_gun ~ timeZ + I(timeZ^2) +
                (1 | intv_area) + (1 | block_id),
              zi ~ 1 + (1 | intv_area) + (1 | block_id)),
           data = d_post,
           family = zero_inflated_negbinomial(),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(normal(0, 0.1), class = b),
                     prior(exponential(2), class = sd),
                     prior(exponential(0.02), class = shape),
                     prior(logistic(0, 1), class = Intercept, dpar = zi),
                     prior(exponential(2), class = sd, dpar = zi)),
           cores = parallel::detectCores(),
           chains = 4, iter = 2000, warmup = 1000, threads = 2,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 15),
           seed = 596692,
           file = "gvff - zinb quad",
           file_refit = "on_change",
           backend = "rstan")

quad = add_criterion(quad, "waic")

lin_random = brm(bf(total_gun ~ timeZ +  
                      (1 + timeZ | intv_area) + (1 + timeZ | block_id),
                    zi ~ 1 + (1 | intv_area) + (1 | block_id)),
                 data = d_post,
                 family = zero_inflated_negbinomial(),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(exponential(0.02), class = shape),
                           prior(lkj(4), class = cor),
                           prior(logistic(0, 1), class = Intercept, dpar = zi),
                           prior(exponential(2), class = sd, dpar = zi)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 181351,
                 file = "gvff - zinb linran",
                 file_refit = "on_change",
                 backend = "rstan")

lin_random = add_criterion(lin_random, "waic")

quad_random = brm(bf(total_gun ~ timeZ + I(timeZ^2) + 
                    (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                    zi ~ 1 + (1 | intv_area) + (1 | block_id)),
                  data = d_post,
                  family = zero_inflated_negbinomial(),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(normal(0, 0.1), class = b),
                            prior(exponential(2), class = sd),
                            prior(exponential(0.02), class = shape),
                            prior(lkj(4), class = cor),
                            prior(logistic(0, 1), class = Intercept, dpar = zi),
                            prior(exponential(2), class = sd, dpar = zi)),
                  cores = parallel::detectCores(),
                  chains = 4, iter = 2000, warmup = 1000, threads = 2,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15),
                  seed = 407013,
                  file = "gvff - zinb quadran",
                  file_refit = "on_change",
                  backend = "rstan")

quad_random = add_criterion(quad_random, "waic")

# Compare (Left column, Table 4)

loo_compare(linear, 
            quad, 
            lin_random,
            quad_random, 
            criterion = "waic")

## Shots fired

linear = brm(bf(cfs_shotsfired ~ timeZ +  
                  (1 | intv_area) + (1 | block_id),
                zi ~ 1 + (1 | intv_area) + (1 | block_id)),
             data = d_post,
             family = zero_inflated_negbinomial(),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(exponential(0.02), class = shape),
                       prior(logistic(0, 1), class = Intercept, dpar = zi),
                       prior(exponential(2), class = sd, dpar = zi)),
             cores = parallel::detectCores(),
             chains = 4, iter = 2000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 551258,
             file = "model fits/sfff - zinb lin",
             file_refit = "on_change",
             backend = "rstan")

linear = add_criterion(linear, "waic")

quad = brm(bf(cfs_shotsfired ~ timeZ + I(timeZ^2) +
                (1 | intv_area) + (1 | block_id),
              zi ~ 1 + (1 | intv_area) + (1 | block_id)),
           data = d_post,
           family = zero_inflated_negbinomial(),
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(normal(0, 0.1), class = b),
                     prior(exponential(2), class = sd),
                     prior(exponential(0.02), class = shape),
                     prior(logistic(0, 1), class = Intercept, dpar = zi),
                     prior(exponential(2), class = sd, dpar = zi)),
           cores = parallel::detectCores(),
           chains = 4, iter = 2000, warmup = 1000, threads = 2,
           control = list(adapt_delta = 0.99,
                          max_treedepth = 15),
           seed = 950067,
           file = "sfff - zinb quad",
           file_refit = "on_change",
           backend = "rstan")

quad = add_criterion(quad, "waic")

lin_random = brm(bf(cfs_shotsfired ~ timeZ +  
                      (1 + timeZ | intv_area) + (1 + timeZ | block_id),
                    zi ~ 1 + (1 | intv_area) + (1 | block_id)),
                 data = d_post,
                 family = zero_inflated_negbinomial(),
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.1), class = b),
                           prior(exponential(2), class = sd),
                           prior(exponential(0.02), class = shape),
                           prior(lkj(4), class = cor),
                           prior(logistic(0, 1), class = Intercept, dpar = zi),
                           prior(exponential(2), class = sd, dpar = zi)),
                 cores = parallel::detectCores(),
                 chains = 4, iter = 2000, warmup = 1000, threads = 2,
                 control = list(adapt_delta = 0.99,
                                max_treedepth = 15),
                 seed = 776398,
                 file = "sfff - zinb linran",
                 file_refit = "on_change",
                 backend = "rstan")

lin_random = add_criterion(lin_random, "waic")

quad_random = brm(bf(cfs_shotsfired ~ timeZ + I(timeZ^2) + 
                       (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                     zi ~ 1 + (1 | intv_area) + (1 | block_id)),
                  data = d_post,
                  family = zero_inflated_negbinomial(),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(normal(0, 0.1), class = b),
                            prior(exponential(2), class = sd),
                            prior(exponential(0.02), class = shape),
                            prior(lkj(4), class = cor),
                            prior(logistic(0, 1), class = Intercept, dpar = zi),
                            prior(exponential(2), class = sd, dpar = zi)),
                  cores = parallel::detectCores(),
                  chains = 4, iter = 2000, warmup = 1000, threads = 2,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15),
                  seed = 574706,
                  file = "sfff - zinb quadran",
                  file_refit = "on_change",
                  backend = "rstan")

quad_random = add_criterion(quad_random, "waic")

# Compare (Right column, Table 4)

loo_compare(linear, 
            quad, 
            lin_random, 
            quad_random, 
            criterion = "waic")

#### Modeling for Inference ####

## Gun Violence

gv_fit = brm(bf(total_gun ~ 
               timeZ + I(timeZ^2) + fs_between + fs_within +
               tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
               (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
               zi ~ 1 + (1 | intv_area) + (1 | block_id)),
             data = d_post,
             family = zero_inflated_negbinomial(),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(exponential(0.02), class = shape),
                       prior(lkj(4), class = cor),
                       prior(logistic(0, 1), class = Intercept, dpar = zi),
                       prior(exponential(2), class = sd, dpar = zi)),
             cores = parallel::detectCores(),
             chains = 4, iter = 4000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 111992,
             file = "zi gun violence fit",
             file_refit = "on_change",
             backend = "rstan")

# Gun Violence Posterior Summaries - Table 5

print(gv_fit, digits = 3) 

## Calls for Service - Shots Fired

sf_fit = brm(bf(cfs_shotsfired ~ 
                  timeZ + I(timeZ^2) + fs_between + fs_within +
                  tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                  (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                zi ~ 1 + (1 | intv_area) + (1 | block_id)),
             data = d_post,
             family = zero_inflated_negbinomial(),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(exponential(0.02), class = shape),
                       prior(lkj(4), class = cor),
                       prior(logistic(0, 1), class = Intercept, dpar = zi),
                       prior(exponential(2), class = sd, dpar = zi)),
             cores = parallel::detectCores(),
             chains = 8, iter = 2500, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 95069,
             file = "zi shots fired fit",
             file_refit = "on_change",
             backend = "rstan")

# Shots Fired Posterior Summaries - Table 5

print(sf_fit, digits = 3)

#### Robustness Checks

### Alternate Dependent Variable

ngv_fit = brm(bf(total_nongun ~ 
                  timeZ + I(timeZ^2) + fs_between + fs_within +
                  tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                  (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                zi ~ 1 + (1 | intv_area) + (1 | block_id)),
             data = d_post,
             family = zero_inflated_negbinomial(),
             prior = c(prior(normal(0, 1), class = Intercept),
                       prior(normal(0, 0.1), class = b),
                       prior(exponential(2), class = sd),
                       prior(exponential(0.02), class = shape),
                       prior(lkj(4), class = cor),
                       prior(logistic(0, 1), class = Intercept, dpar = zi),
                       prior(exponential(2), class = sd, dpar = zi)),
             cores = parallel::detectCores(),
             chains = 4, iter = 4000, warmup = 1000, threads = 2,
             control = list(adapt_delta = 0.99,
                            max_treedepth = 15),
             seed = 593170,
             file = "zi alt dv - gun violence fit",
             file_refit = "on_change",
             backend = "rstan")

# Reverse Time Order

rtogv_fit = brm(bf(rto_total_gun ~ 
                   timeZ + I(timeZ^2) + fs_between + fs_within +
                   tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z +
                   (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                 zi ~ 1 + (1 | intv_area) + (1 | block_id)),
              data = d_post,
              family = zero_inflated_negbinomial(),
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 0.1), class = b),
                        prior(exponential(2), class = sd),
                        prior(exponential(0.02), class = shape),
                        prior(lkj(4), class = cor),
                        prior(logistic(0, 1), class = Intercept, dpar = zi),
                        prior(exponential(2), class = sd, dpar = zi)),
              cores = parallel::detectCores(),
              chains = 4, iter = 4000, warmup = 1000, threads = 2,
              control = list(adapt_delta = 0.99,
                             max_treedepth = 15),
              seed = 337345,
              file = "zi rto - gun violence fit",
              file_refit = "on_change",
              backend = "rstan")

### Compute Differences in Posteriors

library(tidybayes)

# Extract Posteriors

gv_fsw = gv_fit %>%
  spread_draws(b_fs_within) 

ngv_fsw = ngv_fit %>%
  spread_draws(b_fs_within)

rto_gv_fsw = rtogv_fit %>%
  spread_draws(b_fs_within)

## Differences in Posteriors

(gv_fsw$b_fs_within - ngv_fsw$b_fs_within) %>%
  mean_qi(.width = 0.95)

(gv_fsw$b_fs_within - rto_gv_fsw$b_fs_within) %>%
  mean_qi(.width = 0.95)

# Visualize

tibble(
  diff = c(gv_fsw$b_fs_within - ngv_fsw$b_fs_within,
           gv_fsw$b_fs_within - rto_gv_fsw$b_fs_within),
  outcome = rep(c("Alternate Dependent Variable\nNon-Gun Violence", 
                  "Reverse Time Order\nPre-Intervention Gun Violence"), each = 12000)
) %>%
  ggplot(aes(x = diff, fill = after_stat(x < 0))) +
  stat_halfeye(.width = c(0.9, 0.95)) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = expression(atop("Difference in Posterior Distributions", paste("Main Result ", beta[w], " - Comparison ", beta[w]))),
       y = "Posterior Density") +
  facet_wrap(~outcome, nrow = 1) +
  scale_fill_manual(values = c("gray", "deepskyblue"),
                    guide = "none") +
  theme_classic() +
  theme(text = element_text(color = "black", family = "serif"),
        axis.text = element_text(color = "black", size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 12, face = "bold"))

# Figure 3

ggsave("Placebo Posterior Differences.png", width = 6.5, height = 4, units = "in")

#### Stacked and Interaction Version ####

d_post_gv = d_post %>%
  mutate(outcome = "gun violence",
         dv = total_gun)

d_post_ngv = d_post %>%
  mutate(outcome = "non-gun violence",
         dv = total_nongun)

d_post_stack = d_post_gv %>%
  bind_rows(d_post_ngv)

d_post_stack = d_post_stack %>%
  mutate(outcome = factor(outcome, levels = c("non-gun violence", "gun violence")))

stackgv_fit = brm(bf(dv ~ 
                       outcome*(timeZ + I(timeZ^2) + fs_between + fs_within +
                       tx + traffic_sqmiZ + t.popZ + femhhZ + vacantZ + blackZ + hispanicZ + ml1521Z) +
                       (1 + timeZ + I(timeZ^2) | intv_area) + (1 + timeZ + I(timeZ^2) | block_id),
                     zi ~ 1 + outcome + (1 | intv_area) + (1 | block_id)),
                  data = d_post_stack,
                  family = zero_inflated_negbinomial(),
                  prior = c(prior(normal(0, 1), class = Intercept),
                            prior(normal(0, 0.1), class = b),
                            prior(exponential(2), class = sd),
                            prior(exponential(0.02), class = shape),
                            prior(lkj(4), class = cor),
                            prior(logistic(0, 1), class = Intercept, dpar = zi),
                            prior(exponential(2), class = sd, dpar = zi)),
                  cores = parallel::detectCores(),
                  chains = 4, iter = 4000, warmup = 1000, threads = 2,
                  control = list(adapt_delta = 0.99,
                                 max_treedepth = 15),
                  seed = 95197,
                  file = "model fits/stacked dv check",
                  file_refit = "on_change",
                  backend = "rstan")

print(stackgv_fit, digits = 3)
