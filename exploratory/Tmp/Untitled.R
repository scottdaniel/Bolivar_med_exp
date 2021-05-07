load("output/FT_rarefied_w_taxa_glom_df.Rdata")

s24 = ft_w_taxa %>% filter(grepl("Muribaculaceae", Family))

s24_edit = s24 %>% group_by_at(vars(Domain:Species)) %>% summarise(across(starts_with("X"), sum))
    
s24_edit$Group = c(1, 2, 2, 2, 3, 3, 4, 5, 6, 6, 6, 7, 6)

s24_edit2 <- s24_edit %>% group_by_at(vars(Domain:Family, Group)) %>% summarise(across(starts_with("X"), sum))
s24_edit2$Group <- c("unknown", "CAG-873", "metagenome", "Murribaculum", "Porphyromonadaceae bacterium C941", "unclutured", "uncultured Porphromonadaceae")

s24_f <- s24_edit2 %>% pivot_longer(cols = starts_with("X"), names_to = "sample_name", values_to = "Count") %>% filter(Count>0)

s24_f_w <- merge(s24_f, Mapping_work2, by = 0, all.x = T)

ggplot(s24_f_w %>% filter(Body_Site=="Feces"), aes(x = Group, y = Count/10000)) + geom_boxplot() + facet_grid(Year~Urban) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()

ggplot(s24_f_w, aes(x = Group, y = Count/10000)) + geom_boxplot() + facet_grid(Body_Site~Urban) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_log10()

Mapping_work2 %>% group_by(Year, Body_Site, Urban) %>% summarise(N = n())

s24_f_w %>% group_by(Year, Body_Site, Urban) %>% summarise(N = length(unique(Row.names)))
