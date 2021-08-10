library(parallel)
library(MASS)
library(tidyverse)
library(randomForest)

load("output/alpha_imported_all-20190722.Rdata")
load("output/FT_rarefied_full_11142019.Rdata")
load("output/Mapping_MS_Urban.Rdata")

alpha_w_map <- merge(alpha, Mapping_MS_Urban, by = 1)
Work <- alpha_w_map %>% mutate(Urban = ifelse(Ethnicity=="SANEMA", "Low", ifelse(Village=="Fiyakwanha", "Low", ifelse(SampleGroup=="Visitors", "High", "Medium"))) %>% factor(levels = c("Low", "Medium", "High")), Age_num = as.numeric(Age), Age_num = ifelse(is.na(Age_num), as.numeric(gsub(pattern = "_months", replacement = "", Age, ignore.case = T))/12, Age_num), Age_grp1 = ifelse(Age_num>=18, "Adults", "Children"), Age_grp2 = case_when(Age_num<=6 ~ "Age1", Age_num<=12 ~ "Age2", Age_num<=18 ~ "Age3", T ~ "Age4"))

Indices <- colnames(alpha[2:5])
Yrs <- c("2015", "2016")
BodySites <- c("Feces", "Mouth", "Nose", "Right_Arm", "Right_Hand")

run_rf <- function(BB, YY, meta, ft, output){
    Work_dat <- meta %>% filter(Body_Site==BB, Year==YY, SampleGroup=="Villagers")
    Work_biom <- ft[Work_dat$Row.names]
    Work_biom <- Work_biom[rowSums(Work_biom)>0, ]
    rf_x <- t(Work_biom) %>% merge(., Work_dat %>% select(Row.names, Urban, Gender), by.x = 0, by.y = 1) %>% mutate(Urban = as.factor(Urban), Gender = as.factor(Gender)) %>% column_to_rownames(var = "Row.names")
    ts()
    print("Start RF")
    output[[paste(BB, YY, sep = "-")]] <- randomForest(x = rf_x, y = NULL, ntree = 10000, proximity = T, oob.prox = T)
    print("Finish RF")
    ts()
}

rf_fits <- list()
pars <- expand.grid(Yrs, BodySites)
mclapply(seq(10), function(x) run_rf(pars[x, 2], pars[x, 1], Work, biom, rf_fits), mc.cores = 10)

save(rf_fits, "output/rf_fits.Rdata")
