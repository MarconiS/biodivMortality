library(tidyverse)
library(lme4)
library(brms)
library(bayestestR)
library(sf)

sites = c("SERC", "OSBS")
q2 = models = list()
for(site in sites){
data = read_sf(paste("data/",site,"_gap_features.gpkg", sep="")) %>% 
    data.frame %>% select(-one_of("geom")) %>%  filter(dead_count>0) %>%
    filter(raster_value %in% c(41,42,43,52,71, 90,95)) %>% filter(hill_div >0)

fixed_effects = c('surface_dead', 'dead_count', 'dominant_dead', 'dominant_dead_ratio', 'network_length')
random_effects = c('raster_value', 'sci_name')

data[fixed_effects] = scale(data[fixed_effects] )
# 
frml = as.formula(paste("hill_div ~", paste(fixed_effects, collapse = " + "), " + ",
                        paste( "((", paste(fixed_effects, collapse = " + "), ") | ",
                         #"sci_name : raster_value", 
                         random_effects,
                          ")", collapse = " + ")))
#calculate priors assuming that all features but network_length follow an exponential decay distribution
# Prior for 'dead' follows exponential decay
prior = get_prior(frml, data = data, family = lognormal())
prior_surface_dead <- prior(exponential(1), class = 'b', coef = "surface_dead")
prior_dead_count <- prior(exponential(1), class = 'b', coef = "dead_count")
prior_dominant_dead <- prior(exponential(1), class = 'b', coef = "dominant_dead")
prior_dominant_dead_ratio <- prior(exponential(1), class = 'b', coef = "dominant_dead_ratio")

#prior for network_length follows a normal distribution
prior_length = prior(normal(0, 1), class = 'b', coef = "network_length")

# usigng frml as formula, build a brms model
model = brm(frml, data,  family = lognormal(), 
    #prior = c(prior_surface_dead, prior_dead_count, prior_dominant_dead, prior_dominant_dead_ratio, prior_length),
    chains = 2, cores = 2, iter = 6000, 
    threads = threading(6, grainsize = 100), backend = "cmdstanr")


summary(model)
q2[[site]] = describe_posterior(model)
models[[site]] = model
}
saveRDS(q2, "outdir/q2_b_sercosbs_summary.rds")
saveRDS(models, "outdir/models_b_sercosbs_summary.rds")

#describe posterior marginalized by random effects
describe_posterior(model, effects = "all")

# using bayestestr package, summarize rope, significance and direction for each feature
sites = c("HARV", "SERC", "OSBS")
q2 = models = list()
for(site in sites){
data = read_sf(paste("data/",site,"_gap_features.gpkg", sep="")) %>% 
    data.frame %>% select(-one_of("geom")) %>%  filter(dead_count>0) %>%
    filter(raster_value %in% c(41,42,43,52,71, 90,95)) %>% filter(hill_div >0)

fixed_effects = c('surface_dead', 'dead_count', 'dominant_dead', 'dominant_dead_ratio', 'network_length')
fixed_effects = colnames(data)[10:16]
random_effects = c('raster_value', 'sci_name')

data[fixed_effects] = scale(data[fixed_effects] )
# 
frml = as.formula(paste("hill_div ~", paste(fixed_effects, collapse = " + "), " + ",
                        paste( "((", paste(fixed_effects, collapse = " + "), ") | ",
                         #"sci_name : raster_value", 
                         random_effects,
                          ")", collapse = " + ")))
#calculate priors assuming that all features but network_length follow an exponential decay distribution
# Prior for 'dead' follows exponential decay
prior = get_prior(frml, data = data, family = lognormal())
prior_surface_dead <- prior(exponential(1), class = 'b', coef = "surface_dead")
prior_dead_count <- prior(exponential(1), class = 'b', coef = "dead_count")
prior_dominant_dead <- prior(exponential(1), class = 'b', coef = "dominant_dead")
prior_dominant_dead_ratio <- prior(exponential(1), class = 'b', coef = "dominant_dead_ratio")

#prior for network_length follows a normal distribution
prior_length = prior(normal(0, 1), class = 'b', coef = "network_length")

# usigng frml as formula, build a brms model
model = brm(frml, data,  family = lognormal(), 
    #prior = c(prior_surface_dead, prior_dead_count, prior_dominant_dead, prior_dominant_dead_ratio, prior_length),
    chains = 2, cores = 2, iter = 6000, 
    threads = threading(6, grainsize = 100), backend = "cmdstanr")


summary(model)
q2[[site]] = describe_posterior(model)
models[[site]] = model
}
saveRDS(q2, "outdir/q2_b_demogr_summary.rds")
saveRDS(models, "outdir/models_b_demogr_summary.rds")
