library(tidyverse)

load("output/Mapping_MS_Urban.Rdata")
load("output/FT_rarefied_w_taxa_glom.Rdata")
source("ANCOM_updated_code_06_06_2020.R")

Yrs <- c("2015", "2016")
BodySites <- c("Feces", "Mouth", "Nose", "Right_Arm", "Right_Hand")

Mapping_work <- Mapping_MS_Urban %>% mutate(Urban = ifelse(Ethnicity=="SANEMA", "Low", ifelse(Village=="Fiyakwanha", "Low", ifelse(SampleGroup=="Visitors", "High", "Medium"))) %>% factor(levels = c("Low", "Medium", "High")), Age_num = as.numeric(Age), Age_num = ifelse(is.na(Age_num), as.numeric(gsub(pattern = "_months", replacement = "", Age, ignore.case = T))/12, Age_num), Age_grp1 = ifelse(Age_num>=18, "Adults", "Children"))
Mapping_work$Age_num[Mapping_work$Age=="4_Days"] <- 0
Mapping_work$Age_num[Mapping_work$Age=="1_day"] <- 0
Mapping_work$Age_grp1[Mapping_work$Age_num==0] <- "Children"

Mapping_work2 <- Mapping_work %>% mutate(Age_grp2 = case_when(Age_num<=3 ~ "Age_0-3", Age_num<=8 ~ "Age_3-8", Age_num<=18 ~ "Age_8-18", T ~ "Adults") %>% factor(levels = c("Age_0-3", "Age_3-8", "Age_8-18", "Adults")))

ancom_wrap <- function(featuretable, metadata, ...){
    ft <- featuretable # feature table is in a format that a column of taxa with additional columns of counts from each samples
    mt <- metadata

    # filter feature table to remove features that appeared in less than 10% of samples
    sample_names <- colnames(ft)[-1]
    n_samples <- length(sample_names)
    ft_flt <- ft[rowSums(ft[, -1]>0)>=0.1*n_samples, ]

    ## create ancom format feature table
    ft_ancom <- ft_flt[, -1] %>% t() %>% as.data.frame() %>% rownames_to_column(var = "Sample.ID")
    colnames(ft_ancom) <- c("Sample.ID", ft_flt$Taxa) %>% make.names()


    # filter metadata
    mt_ancom <- mt %>% filter(SampleID %in% ft_ancom$Sample.ID) %>% rename(Sample.ID = SampleID)

    ancom_mod <- ANCOM.main(OTUdat = ft_ancom, Vardat = mt_ancom, ...)
    return(ancom_mod)

}

# running ancom this process can be quite long, each process taking about 20 minutes
Ancom_mf_lst_2 <- Mapping_work2 %>% filter(SampleGroup=="Villagers", Body_Site %in% BodySites) %>% group_by(Year, Body_Site) %>% group_split()

Ancom_rsl_2 <- parallel::mclapply(Ancom_mf_lst_2, function(x){
    list(x[1, c("Year", "Body_Site")], ancom_wrap(ft_glom$Genus %>% select(Taxa, any_of(x$SampleID)), x, adjusted = T, main.var = "Urban", adj.formula = "Gender+Age_num", repeat.var = NULL, longitudinal = F, random.formula = NULL, multcorr = 2, sig = 0.05, prev.cut = 1))
    }, mc.cores = 20)

save(Ancom_rsl_2, file = "output/Ancom_rsl_g_YYBB.Rdata")

