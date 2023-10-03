
sites = c("HARV", "SERC", "OSBS")
q1 = dtt = list()
for(site in sites){
    data_no_dead = read_sf(paste("data/",site,"_no_dead_gap_features.gpkg", sep="")) %>% 
        data.frame %>% select(-one_of("geom")) %>%  filter(dead_count==0) %>%
            filter(raster_value %in% c(41,42,43,52,71, 90,95)) %>% filter(hill_div >0) %>%
        mutate(deads_alive = "alive") %>% 
        select(hill_div, sci_name, raster_value, deads_alive)

    data_dead = read_sf(paste("data/",site,"_gap_features.gpkg", sep="")) %>% 
        data.frame %>% select(-one_of("geom")) %>%  filter(dead_count>0) %>%
            filter(raster_value %in% c(41,42,43,52,71, 90,95)) %>% filter(hill_div >0) %>%
        mutate(deads_alive = "dead") %>% 
        select(hill_div, sci_name, raster_value, deads_alive)

    # from data_no_dead, and data_dead select hill, sci_name and raster_value
    data = rbind.data.frame(data_no_dead, data_dead) 
    data$raster_value = as.character(data$raster_value)
    data$deads_alive = as.factor(data$deads_alive)
    data =  data %>% filter(raster_value %in% c(41, 42, 43, 52, 90, 95))
    # for each raster_value and sci_name, test whether there is a significant difference between deads_alive and hill_div
    list_p_values = list()
    ii = 0
    for(nlcd in unique(data$raster_value)){
        for(sname in unique(data$sci_name)){
            ii = ii + 1
            data_subset = data %>% dplyr::filter(sci_name == sname, raster_value == nlcd)
            # if less than 3 dead or less than 3 alive, skip
            if(length(data_subset[["deads_alive"]][data_subset$deads_alive=="dead"]) < 3){
                next
            }
            if(length(data_subset[["deads_alive"]][data_subset$deads_alive=="alive"]) < 3){
                next
            }
            # resample dataset to be balanced across deads_alive factor values, for each sci_name
            # Determine the minimum count between "dead" and "alive"
            min_count <- min(
                sum(data_subset$deads_alive == "dead"),
                sum(data_subset$deads_alive == "alive")
            )

            # Now sample min_count instances from each group to balance the dataset
            balanced_dataset <- data_subset %>%
                group_by(deads_alive) %>%
                sample_n(min_count) %>%
                ungroup() 


            # Now you can proceed with your original code...
            data_subset = balanced_dataset %>% group_by(deads_alive)
            list_p_values[[ii]]  = test_if_dead(balanced_dataset)        
        }
    }
    list_p_values = do.call(rbind, list_p_values)
    q1[[site]] = data.frame(list_p_values, site = site, mean_hill_div = mean(data$hill_div))
        dtt[[site]] = data.frame(data, site = site)
}
q1_results = do.call(rbind, q1)
q1_results
dtt = do.call(rbind, dtt)

# plot the distribution of hill for  each site, color coded by deads_alive
ggplot(dtt, 
aes(y=hill_div, fill=deads_alive, x = deads_alive)) + 
geom_violin(alpha=0.2) + 
facet_wrap(~site, ncol=3) +
ggsci::scale_fill_jco() +
theme_bw() 
ggsave("figures/hill_dist.png", width=10, height=3, units="in", dpi=300)

q1_results["raster_value"] = ifelse(q1_results$raster_value == "41", "Deciduous Forest", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "42", "Evergreen Forest", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "43", "Mixed Forest", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "52", "Shrubs", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "90", "Woody Wetlands", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "95", "Emergent Herbaceous Wetlands", q1_results$raster_value)

# scatterplot of r2 values for each site and each raster_value type
q1_results %>% filter(raster_value %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest",  "Woody Wetlands")) %>%
ggplot(aes(y=raster_value, x=r2, color = site)) +
geom_point() + geom_jitter() + ggsci::scale_color_jco() + theme_bw() + theme(legend.position = "bottom") 
ggsave("figures/q1_r2.png", width=6, height=3, units="in", dpi=300)


q1_results %>% filter(raster_value %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest",  "Woody Wetlands")) %>%
ggplot(aes(y=raster_value, x=slope, color = site)) + geom_vline(xintercept = 0., linetype="dashed") +
geom_point() + geom_jitter() + ggsci::scale_color_jco() + theme_bw() + theme(legend.position = "bottom") 
ggsave("figures/q1_slopes.png", width=6, height=3, units="in", dpi=300)

#barplot of the raster_values with NA vs not NA
q1_results %>%
    filter(raster_value %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest",  "Woody Wetlands")) %>%
    group_by(raster_value, site) %>%
    # summarize the ratio of NA to not NA in r2
    summarise(n = sum(!is.na(r2))/n()) %>%
    ggplot(aes(y=n, x=raster_value, fill = site)) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), width=0.8) +
    ggsci::scale_fill_jco() + ylim(0,1)+
    theme_bw() +
    theme(legend.position = "bottom")
ggsave("figures/ratio_significant.png", width=6, height=3, units="in", dpi=300)

q1_results %>%
    filter(raster_value %in% c("Deciduous Forest"), site == "OSBS")

q1_results %>%
    filter(raster_value %in% c("Deciduous Forest", "Evergreen Forest", "Mixed Forest",  "Woody Wetlands")) %>%
    group_by(raster_value, sci_name) %>%
    # summarize the ratio of NA to not NA in r2
    summarise(n = sum(!is.na(r2))/n()) %>% filter(n>0) %>%
    ggplot(aes(x=n, y=sci_name)) + geom_point()+
    theme_bw() +
    theme(legend.position = "bottom")
ggsave("figures/ratio_significant_sci_name.png", width=6, height=5, units="in", dpi=300)
