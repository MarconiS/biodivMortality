library(tidyverse)
library(lme4)
library(sf)
library(ggplot2)

test_if_dead = function(x){
    # x is a dataframe with columns hill_div, sci_name, raster_value, deads_alive
    # test if there is a significant difference between deads_alive and hill_div
    # return a dataframe with columns sci_name, raster_value, p_value
    # p_value is the p_value of the t-test
    # sci_name and raster_value are the same for all rows
    # p_value is the same for all rows
    # if there is no significant difference, return NA
    # if there is a significant difference, return the p_value
    # if there is only one row, return NA
    # if there are no rows, return NA
    if(nrow(x) == 0){
        return(data.frame(sci_name = NA, raster_value = NA, p_value = NA))
    }
    if(nrow(x) == 1){
        return(data.frame(sci_name = NA, raster_value = NA, p_value = NA))
    }
    if(nrow(x) == 2){
        return(data.frame(sci_name = NA, raster_value = NA, p_value = NA))
    }
    if(nrow(x) > 2){
        p_value = summary(lm(x$hill_div ~ x$deads_alive))
        f = p_value$fstatistic        
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        if(p > 0.05){
            return(data.frame(sci_name = x$sci_name[1], raster_value = x$raster_value[1], p_value = NA, r2 = NA, intercept = NA, slope = NA))
        }else{
            return(data.frame(sci_name = x$sci_name[1], raster_value = x$raster_value[1], p_value = p, r2 = p_value$r.squared, intercept = p_value$coefficients[1], slope = p_value$coefficients[2]))
        }
    }
}


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
            list_p_values[[ii]]  = test_if_dead(data_subset)        
        }
    }
    list_p_values = do.call(rbind, list_p_values)
    q1[[site]] = data.frame(list_p_values, site = site)
        dtt[[site]] = data.frame(data, site = site)
}
q1_results = do.call(rbind, q1)
dtt = do.call(rbind, dtt)

q1_results["significant"] = ifelse(!is.na(q1_results$p_value), "yes", "no")
q1_results["sign"] = ifelse(q1_results$slope>0, "positive", "negative")


#ggplot barplot for each raster_value, with x being the proportion of significant, grouped by site
q1_results %>% group_by(raster_value, site, significant) %>% summarize(n = n()) %>% 
    ggplot(aes(x = site, y = n, fill = significant)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(.~raster_value) + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "NLCD code", y = "Number of species") + 
    ggtitle("Number of species with significant difference between dead and alive") + 
    ggsci::scale_fill_jco() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./figures/q1_significant_mortality_effect.png", width = 10, height = 10, units = "in")
#ggplot to show the proportion of positive and negative slopes, grouped by site and NLCD code
q1_results %>% group_by(raster_value, site, sign) %>% summarize(n = n()) %>% 
    ggplot(aes(x = raster_value, y = n, fill = sign)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(.~site) + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "NLCD code", y = "Number of species") + 
    ggtitle("Number of species with positive and negative slopes") + 
    ggsci::scale_fill_jco() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./figures/q1_direction_mortality_effect.png", width = 10, height = 3, units = "in")

# plot the number of times a species has a significant slope, faceted per site
q1_results %>% group_by(sci_name, site, significant) %>% summarize(n = sum(!is.na(p_value))) %>% 
    ggplot(aes(x = sci_name, y = n)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(.~site, nrow=3) + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "Site", y = "Number of significant slopes") + 
    ggtitle("Number of significant slopes per species") + 
    ggsci::scale_fill_jco() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./figures/q1_species_significance_effect.png", width = 10, height = 10, units = "in")

# plot the number of times a species has a significant slope, faceted per site
q1_results %>% group_by(sci_name, site, sign) %>% summarize(n = sum(!is.na(p_value))) %>% 
    ggplot(aes(x = site, y = n, fill = sign)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(.~sci_name) + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "Site", y = "Number of significant slopes") + 
    ggtitle("Number of significant slopes per species") + 
    ggsci::scale_fill_jco() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("./figures/q1_species_direction_effect.png", width = 10, height = 10, units = "in")

# from raster_value change 41 in Deciduous Forest to 42 in Evergreen Forest, and 43 in Mixed Forest, 52 Shrubs, 90 woody wetlands, 95 Emergent Herbaceous Wetlands
q1_results["raster_value"] = ifelse(q1_results$raster_value == "41", "Deciduous Forest", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "42", "Evergreen Forest", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "43", "Mixed Forest", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "52", "Shrubs", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "90", "Woody Wetlands", q1_results$raster_value)
q1_results["raster_value"] = ifelse(q1_results$raster_value == "95", "Emergent Herbaceous Wetlands", q1_results$raster_value)

q1_results %>% group_by(raster_value, site, significant) %>% summarize(n = n()) %>% 
    ggplot(aes(x = site, y = n, fill = significant)) + geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(.~raster_value) + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "NLCD code", y = "Number of species") + 
    ggtitle("Number of species with significant difference between dead and alive") + 
    ggsci::scale_fill_jco() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#plot barplot with: x being raster_value + site, y being the relative proportion of significant slopes, on top of each bar, geom_text the number of significant slopes/total per group
q1_results %>% group_by(raster_value, site, significant) %>% summarize(n = n()) %>% 
    ungroup() %>% group_by(raster_value, site) %>% mutate(perc = n/sum(n), tot = sum(n)) %>%
    filter(significant == "yes") %>%
    ggplot(aes(x = paste(raster_value, site, sep = "\n"), y = perc)) + 
    # on top of column print paste(perc, "/", tot, sep = "")
    geom_col() + geom_text(aes(label = paste(n, "/", tot, sep = "")), hjust = -1.5) +
    geom_bar(stat = "identity") + 
    theme_bw() + theme(legend.position = "bottom") +
    labs(x = "NLCD code", y = "Proportion of species") + 
    ggtitle("Number of species with significant difference between dead and alive") + 
    ggsci::scale_fill_jco() + coord_flip() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 20))
ggsave("./figures/q1_proportion_of_species.png", width = 10, height = 10, units = "in")


# plot the distribution of hill diversity for each site, faceted by raster_value, colored by deads_alive
dtt %>%      
    filter(raster_value %in% c(41,42,43,52,71, 90,95)) %>% filter(hill_div >0) %>% 
    ggplot(aes(x = hill_div, fill = deads_alive)) + 
    geom_density() + 
    facet_wrap(site~., scales = "free") + theme_bw() + theme(legend.position = "bottom") +
    labs(x = "Site", y = "Number of species") + 
    ggtitle("Number of species with significant difference between dead and alive") + 
    ggsci::scale_fill_jco() + 
    # increase size of fonts in plot
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 20)) 

ggsave("./figures/dist_per_site.png", width = 10, height = 3, units = "in")
