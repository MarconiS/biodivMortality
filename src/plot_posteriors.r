library(bayestestR)
odels = readRDS("outdir/models_b_sercosbs_summary.rds")
hrv = readRDS("outdir/models_b_summary.rds")
models[["HARV"]] = hrv[["HARV"]]
q2_and_3 = list()
for(site in sites){
# get posterior description for each feature
q2_and_3[[site]] = data.frame(site=site, describe_posterior(models[[site]], effects = "all"))
# append site name to each feature
}
results_posterior = bind_rows(q2_and_3)


# load library
library(ggplot2)

# iterate over all fixed effect parameters to check if the distributions have any unexpected behavior
for (param in colnames(posterior_samples)) {
  if (!grepl("tercept", param)) {  # fixed effects in brms start with "b_"
    # plot posterior
    ggplot(posterior_samples, aes_string(param)) +
      geom_histogram(aes(y=..density..), bins=30, fill='skyblue', alpha=0.7) +
      geom_density(color="red", size=1.2) +
      theme_classic() +
      labs(x=param, y="Density", title=paste("Posterior distribution of", param)) +
      theme(plot.title = element_text(hjust = 0.5))
    ggsave(paste0("figures/posteriors/",param, ".png"), width=6, height=4, dpi=300)
    # print plot immediately in a loop
    print(last_plot())
  }
}

# plot the CI for each feature starting with b_ and whose Effects is "fixed", colored by site

# Create jitter
set.seed(123)
jitter_val <- ggplot2::position_jitter(width = 0.2, seed = 123)
# factor site so that they appear in the follwong order c("OSBS, "SERC", "HARV")
results_posterior$site <- factor(results_posterior$site, levels = c("OSBS", "SERC", "HARV"))
results_posterior %>% 
  filter(grepl("b_", Parameter), !grepl("Intercept", Parameter), Effects == "fixed") %>%
  mutate(grp =  interaction(site, Parameter)) %>%
  ggplot(aes(x=CI_low, xend=CI_high, y=grp, yend=grp)) +
  geom_vline(xintercept=0, color="grey", size=1.5) +
  geom_point(aes(x=Median, color=site), size=3, alpha=0.5) +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high, color = site), height=0.1) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Estimate", y="Site", title="Posterior for each feature", 
       subtitle="ROPE = Region of Practical Equivalence")
ggsave("figures/landscape_level_effects.png", width=6, height=4, dpi=300)

for(tx in unique(data$sci_name)){
  # substitute spaces in taxa with points
  taxa = gsub(" ", ".", tx)
  results_posterior %>% 
  filter(Parameter %like% taxa) %>%
  #mutate Parameter so that it only contains the content after the last comma
  mutate(Parameter = gsub(".*\\,", "", Parameter)) %>%
  mutate(Parameter = substr(Parameter, 1, nchar(Parameter)-1)) %>%
  mutate(grp =  interaction(site, Parameter)) %>%
  ggplot(aes(x=CI_low, xend=CI_high, y=grp, yend=grp)) +
  geom_vline(xintercept=0, color="grey", size=1.5) +
  geom_point(aes(x=Median, color=site), size=3, alpha=0.5) +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high, color = site), height=0.1) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Estimate", y="Site", title="Posterior for each feature", 
       subtitle=taxa)
  ggsave(paste0("figures/patch/",tx, "_level_effects.png"), width=6, height=4, dpi=300)

}

for(tx in unique(data$raster_value)){
  # substitute spaces in taxa with points
  taxa = gsub(" ", ".", tx)
  results_posterior %>% 
  filter(Parameter %like% taxa) %>%
  #mutate Parameter so that it only contains the content after the last comma
  mutate(Parameter = gsub(".*\\,", "", Parameter)) %>%
  mutate(Parameter = substr(Parameter, 1, nchar(Parameter)-1)) %>%
  mutate(grp =  interaction(site, Parameter)) %>%
  ggplot(aes(x=CI_low, xend=CI_high, y=grp, yend=grp)) +
  geom_vline(xintercept=0, color="grey", size=1.5) +
  geom_point(aes(x=Median, color=site), size=3, alpha=0.5) +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high, color = site), height=0.1) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Estimate", y="Site", title="Posterior for each feature", 
       subtitle=taxa)
  ggsave(paste0("figures/ecosystem_",tx, "_level_effects.png"), width=6, height=4, dpi=300)

}

#create raster_value column extracting the number from the Parameter column
res_ecosys = results_posterior
results_posterior$raster_value = gsub("[^0-9]", "", res_ecosys$Parameter)
res_ecosys = res_ecosys %>% 
    filter(raster_value %in% c("41", "42", "43", "52", "90", "95"))
res_ecosys["raster_value"] = ifelse(res_ecosys$raster_value == "41", "Deciduous Forest", res_ecosys$raster_value)
res_ecosys["raster_value"] = ifelse(res_ecosys$raster_value == "42", "Evergreen Forest", res_ecosys$raster_value)
res_ecosys["raster_value"] = ifelse(res_ecosys$raster_value == "43", "Mixed Forest", res_ecosys$raster_value)
res_ecosys["raster_value"] = ifelse(res_ecosys$raster_value == "52", "Shrubs", res_ecosys$raster_value)
res_ecosys["raster_value"] = ifelse(res_ecosys$raster_value == "90", "Woody Wetlands", res_ecosys$raster_value)
res_ecosys["raster_value"] = ifelse(res_ecosys$raster_value == "95", "Emergent Herbaceous Wetlands", res_ecosys$raster_value)

# get barplot of how many ecosystems are significant for each feature and site
res_ecosys %>% 
filter(!Parameter %like% "ntercept") %>%
  group_by(site, raster_value) %>%
  #mutate Parameter so that it only contains the content after the last comma
  mutate(Parameter = gsub(".*\\,", "", Parameter)) %>%
  mutate(Parameter = substr(Parameter, 1, nchar(Parameter)-1)) %>%
  mutate(grp =  interaction(site, Parameter)) %>%
  ggplot(aes(x=CI_low, xend=CI_high, y=grp, yend=grp)) +
  geom_vline(xintercept=0, color="grey", size=1.5) +
  geom_point(aes(x=Median, color=raster_value), size=3, alpha=0.5, position = jitter_val) +
  geom_errorbarh(aes(xmin=CI_low, xmax=CI_high, color = raster_value), height=0.1, position = jitter_val) +
  ggsci::scale_color_jco() +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Estimate", y="Site", title="Posterior for each feature", 
       subtitle=taxa)
  ggsave(paste0("figures/ecosystem_",tx, "_level_recruit_effects.png"), width=6, height=4, dpi=300)

}


res_ecosys = res_ecosys %>% mutate(signifi = ifelse(pd < 0.95, "no", "yes"))
res_ecosys[["rs_site"]] = paste0(res_ecosys$site, "_", res_ecosys$raster_value)
library(ggplot2)

res_ecosys %>%
    filter(raster_value %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest",  "Woody Wetlands")) %>%
    group_by(raster_value, site) %>%
    # summarize the ratio of signifi=="yes" to total
    summarise(n = sum(signifi=="yes")/n()) %>%
    ggplot(aes(y=raster_value, x=site, fill = n)) + 
    geom_raster() +
    geom_text(data = . %>% filter(n > 0), aes(label=sprintf("%.2f", n)), vjust=1.5, size=3) +
    scale_fill_gradient(low = "white", high = "blue", na.value = "white") +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(paste0("figures/only_looking_at_dead_trees/ecosys_significance.png"), width=3, height=4, dpi=300)
