#################################
### JERONIMO & VARASSIN, 2023 ###
#################################

#### Work directory ####
setwd("X:/RDirectory/Masters")
getwd()



#### Load packages ####
if(!require('DHARMa')) install.packages('DHARMa')
library('DHARMa')

if(!require('ggcorrplot')) install.packages('ggcorrplot')
library('ggcorrplot')

if(!require('MuMIn')) install.packages('MuMIn')
library('MuMIn')

if(!require('rstatix')) install.packages('rstatix')
library('rstatix')

if(!require('tidyverse')) install.packages('tidyverse') # data wrangling
library('tidyverse')

if(!require('vegan')) install.packages('vegan') # for Mantel test
library('vegan')

if(!require('visreg')) install.packages('visreg') 
library('visreg')



#### Load data frames ####
# Farm level data
data.prop <- read.csv("Data/Data_Local.csv")

data.prop <- data.prop %>% 
  mutate_if(is.character, as.factor) %>% 
  mutate(ENNF_MN_1KM = as.numeric(ENNF_MN_1KM))

str(data.prop)

View(data.prop)

# Municipal level data
data.muni <- read.csv("Data/Data_Regional.csv")

data.muni <- data.muni %>% 
  mutate_if(is.character, as.factor) 

str(data.muni)

View(data.muni)



#### Mantel test for spatial self-correlation (Farm level only) ####
# Coordinates matrix
data.pcoo <- data.prop %>%
  select(Latitude, Longitude)

str(data.pcoo)

pcoo.std <- decostand(data.pcoo, method="standardize")

pcoo.ecl <- vegdist(pcoo.std, method="euclidean")

str(pcoo.ecl)

# Production (kg) matrix
data.ppkg <- data.prop %>%
  select(Prod_kg)

str(data.ppkg)

ppkg.std <- decostand(data.ppkg, method="standardize")

ppkg.ecl <- vegdist(ppkg.std, method="euclidean")

str(ppkg.ecl)

# Landscape metrics matrix
data.pmtr <- data.prop %>%
  select(PLANDF_5KM, PDF_5KM, PROXF_MN_5KM, ENNF_MN_5KM, CONNECTF_5KM, SHEI_5KM, II_5KM)

str(data.pmtr)

pmtr.std <- decostand(data.pmtr, method="standardize")

pmtr.ecl <- vegdist(pmtr.std, method="euclidean")

str(pmtr.ecl)

# Mantel test production
test.MTL01 <- mantel(pcoo.ecl, ppkg.ecl, method="pearson", permutations=9999)

test.MTL01

data.MTL01 <- cbind(pcoo.ecl, ppkg.ecl)

str(data.MTL01)

write.table(data.MTL01, "Data/data.MTL01.csv")

# Mantel test landscape metrics
test.MTL02 <- mantel(pcoo.ecl, pmtr.ecl, method="pearson", permutations=9999)

test.MTL02

# Mantel plot landscape metrics
data.MTL02 <- cbind(pcoo.ecl, pmtr.ecl)

write.table(data.MTL02, "Data/data.MTL02.csv")



# Pearson's correlation test of predictor variables

## Farm Level
data.COR01 <- data.prop %>% 
  select(PLANDF_5KM, PDF_5KM, PROXF_MN_5KM, ENNF_MN_5KM, CONNECTF_5KM, SHEI_5KM, II_5KM) %>% 
  rename('F. Cover Percentage'=PLANDF_5KM,
         'F. Patches Density'=PDF_5KM,
         'F. Proximity'=PROXF_MN_5KM,
         'F. Isolation'=ENNF_MN_5KM,
         'F. Connectance'=CONNECTF_5KM,
         'L. Diversity'=SHEI_5KM,
         'L. Intensity'=II_5KM)

test.COR01 <- cor(data.COR01)

p.mat <- cor_pmat(data.COR01)

test.COR01



#### Shapiro-Wilk test for normality ####
data.prop %>%
  shapiro_test(ProdArea_kg)



#### Variable selection #### 
#### Farm level effect raidus 
### FOREST PERCENTAGE OF LANDSCAPE 
v1.1KM <- lm(data=data.prop, ProdArea_kg~PLANDF_1KM)
summary(v1.1KM) #R?=0.1324
v1.2KM <- lm(data=data.prop, ProdArea_kg~PLANDF_2KM)
summary(v1.2KM) #R?=0.1388
v1.3KM <- lm(data=data.prop, ProdArea_kg~PLANDF_3KM)
summary(v1.3KM) #R?=0.1435 ***
v1.4KM <- lm(data=data.prop, ProdArea_kg~PLANDF_4KM)
summary(v1.4KM) #R?=0.1189
v1.5KM <- lm(data=data.prop, ProdArea_kg~PLANDF_5KM)
summary(v1.5KM) #R?=0.1038

### FOREST PATCH DENSITY
v2.1KM <- lm(data=data.prop, ProdArea_kg~PDF_1KM)
summary(v2.1KM) #R?=0.1215
v2.2KM <- lm(data=data.prop, ProdArea_kg~PDF_2KM)
summary(v2.2KM) #R?=0.2725 ***
v2.3KM <- lm(data=data.prop, ProdArea_kg~PDF_3KM)
summary(v2.3KM) #R?=0.2348
v2.4KM <- lm(data=data.prop, ProdArea_kg~PDF_4KM)
summary(v2.4KM) #R?=0.1701
v2.5KM <- lm(data=data.prop, ProdArea_kg~PDF_5KM)
summary(v2.5KM) #R?=0.1880

# FOREST PROXIMITY INDEX
v3.1KM <- lm(data=data.prop, ProdArea_kg~PROXF_MN_1KM)
summary(v3.1KM) #R?=0.0489 ***
v3.2KM <- lm(data=data.prop, ProdArea_kg~PROXF_MN_2KM)
summary(v3.2KM) #R?=0.0186
v3.3KM <- lm(data=data.prop, ProdArea_kg~PROXF_MN_3KM)
summary(v3.3KM) #R?=0.0062
v3.4KM <- lm(data=data.prop, ProdArea_kg~PROXF_MN_4KM)
summary(v3.4KM) #R?=0.0120
v3.5KM <- lm(data=data.prop, ProdArea_kg~PROXF_MN_5KM)
summary(v3.5KM) #R?=0.0168

# FOREST EUCLIDEAN NEXT NEIGHBOR
v4.1KM <- lm(data=data.prop, ProdArea_kg~ENNF_MN_1KM)
summary(v4.1KM) #R?=0.2257
v4.2KM <- lm(data=data.prop, ProdArea_kg~ENNF_MN_2KM)
summary(v4.2KM) #R?=0.3565
v4.3KM <- lm(data=data.prop, ProdArea_kg~ENNF_MN_3KM)
summary(v4.3KM) #R?=0.3233
v4.4KM <- lm(data=data.prop, ProdArea_kg~ENNF_MN_4KM)
summary(v4.4KM) #R?=0.3963 ***
v4.5KM <- lm(data=data.prop, ProdArea_kg~ENNF_MN_5KM)
summary(v4.5KM) #R?=0.3918

# FOREST CONNECTIVITY INDEX 
v5.1KM <- lm(data=data.prop, ProdArea_kg~CONNECTF_1KM)
summary(v5.1KM) #R?=0.1732 ***
v5.2KM <- lm(data=data.prop, ProdArea_kg~CONNECTF_2KM)
summary(v5.2KM) #R?=0.0089
v5.3KM <- lm(data=data.prop, ProdArea_kg~CONNECTF_3KM)
summary(v5.3KM) #R?=0.1083
v5.4KM <- lm(data=data.prop, ProdArea_kg~CONNECTF_4KM)
summary(v5.4KM) #R?=0.1427 
v5.5KM <- lm(data=data.prop, ProdArea_kg~CONNECTF_5KM)
summary(v5.5KM) #R?=0.04598

# LANDSCAPE SHANNON EVENNESS INDEX
v6.1KM <- lm(data=data.prop, ProdArea_kg~SHEI_1KM)
summary(v6.1KM) #R?=0.2033
v6.2KM <- lm(data=data.prop, ProdArea_kg~SHEI_2KM)
summary(v6.2KM) #R?=0.2311 ***
v6.3KM <- lm(data=data.prop, ProdArea_kg~SHEI_3KM)
summary(v6.3KM) #R?=0.2125
v6.4KM <- lm(data=data.prop, ProdArea_kg~SHEI_4KM)
summary(v6.4KM) #R?=0.2126
v6.5KM <- lm(data=data.prop, ProdArea_kg~SHEI_5KM)
summary(v6.5KM) #R?=0.1980

# LANDSCAPE INTENSITY INDEX
v7.1KM <- lm(data=data.prop, ProdArea_kg~II_1KM)
summary(v7.1KM) #R?=0.3877
v7.2KM <- lm(data=data.prop, ProdArea_kg~II_2KM)
summary(v7.2KM) #R?=0.3996 ***
v7.3KM <- lm(data=data.prop, ProdArea_kg~II_3KM)
summary(v7.3KM) #R?=0.3915
v7.4KM <- lm(data=data.prop, ProdArea_kg~II_4KM)
summary(v7.4KM) #R?=0.3567
v7.5KM <- lm(data=data.prop, ProdArea_kg~II_5KM)
summary(v7.5KM) #R?=0.3538



#### GLM MODEL SELECTION FARM LEVEL ####
# Analysis data frame
data.GLM01 <- data.prop %>% 
  select(ProdArea_kg, PLANDF_3KM, PDF_2KM, ENNF_MN_4KM, PROXF_MN_1KM, CONNECTF_4KM,
         SHEI_2KM, II_2KM) %>% 
  rename('PROD'=ProdArea_kg,
         'CVRG'=PLANDF_3KM,
         'DNST'=PDF_2KM,
         'ISLT'=ENNF_MN_4KM,
         'PRXM'=PROXF_MN_1KM,
         'CNNC'=CONNECTF_4KM, 
         'SHEI'=SHEI_2KM,
         'IIND'=II_2KM)

# Models
# NULL MODEL
M0 <- glm(data=data.GLM01, formula=PROD~1, family="gaussian")
summary(M0)

# UNIV MODELS
M1 <- glm(data=data.GLM01, formula=PROD~CVRG, family="gaussian")
summary(M1)
M2 <- glm(data=data.GLM01, formula=PROD~DNST, family="gaussian")
summary(M2)
M3 <- glm(data=data.GLM01, formula=PROD~ISLT, family="gaussian")
summary(M3)
M4 <- glm(data=data.GLM01, formula=PROD~PRXM, family="gaussian")
summary(M4)
M5 <- glm(data=data.GLM01, formula=PROD~CNNC, family="gaussian")
summary(M5)
M6 <- glm(data=data.GLM01, formula=PROD~SHEI, family="gaussian")
summary(M6)
M7 <- glm(data=data.GLM01, formula=PROD~IIND, family="gaussian")
summary(M7)



# Univ models selection
UVAR.SLCT <- model.sel(M0, M1, M2, M5, M6, M7)

UVAR.SLCT

write.table(UVAR.SLCT, "Data/UnivModels.csv")

# Best fitted model evaluation
testResiduals(M7)

testResiduals(M2)

testResiduals(M6)

# MULTIV MODELS
MA <- glm(data=data.GLM01, formula=PROD~IIND*CVRG, family="gaussian")
summary(MA)

MB <- glm(data=data.GLM01, formula=PROD~IIND*DNST, family="gaussian")
summary(MB)

MC <- glm(data=data.GLM01, formula=PROD~IIND*PRXM, family="gaussian")
summary(MC)

MD <- glm(data=data.GLM01, formula=PROD~IIND*CNNC, family="gaussian")
summary(MD)

ME <- glm(data=data.GLM01, formula=PROD~IIND*SHEI, family="gaussian")
summary(ME)
         
# Univ models selection
MVAR.SLCT <- model.sel(M0, MA, MB, MD, ME)

MVAR.SLCT

write.table(MVAR.SLCT, "Data/MulvModels.csv")

# Best fitted model evaluation
testResiduals(MB)

testResiduals(MD)

testResiduals(ME)



############################
##### Municipality level ###
############################
data.COR02 <- data.muni %>% 
  select(PLANDF, PDF, PROXF_MN, ENNF_MN, CONNECTF, SHEI, II) %>% 
  rename('F. Cover Percentage'=PLANDF,
         'F. Patches Density'=PDF,
         'F. Proximity'=PROXF_MN,
         'F. Isolation'=ENNF_MN,
         'F. Connectance'=CONNECTF,
         'L. Diversity'=SHEI,
         'L. Intensity'=II)

test.COR02 <- cor(data.COR02)

p.mat2 <- cor_pmat(data.COR02)

test.COR02



# Shapiro-Wilk test for normality
data.muni %>%
  shapiro_test(ProdArea_kg)

ggplot(data.muni, aes(x=ProdArea_kg))+
  geom_histogram(fill="#91091E",colour="black", bins=10)

### Data set GLM City Level
data.GLM02 <- data.muni %>% 
  select(ProdArea_kg, PLANDF, PDF, ENNF_MN, PROXF_MN, CONNECTF, SHEI, SHDI, II) %>% 
  rename('PROD'=ProdArea_kg,
         'CVRG'=PLANDF,
         'DNST'=PDF,
         'ISLT'=ENNF_MN,
         'PRXM'=PROXF_MN,
         'CNNC'=CONNECTF, 
         'SHEI'=SHEI,
         'SHDI'=SHDI,
         'IIND'=II)

# UNIV MODELS
# NULL MODEL
C0 <- glm(data=data.GLM02, formula=PROD~1, family="gaussian")
summary(C0)
C1 <- glm(data=data.GLM02, formula=PROD~CVRG, family="gaussian")
summary(C1)
C2 <- glm(data=data.GLM02, formula=PROD~DNST, family="gaussian")
summary(C2)
C3 <- glm(data=data.GLM02, formula=PROD~ISLT, family="gaussian")
summary(C3)
C4 <- glm(data=data.GLM02, formula=PROD~PRXM, family="gaussian")
summary(C4)
C5 <- glm(data=data.GLM02, formula=PROD~CNNC, family="gaussian")
summary(C5)
C6 <- glm(data=data.GLM02, formula=PROD~SHEI, family="gaussian")
summary(C6)
C7 <- glm(data=data.GLM02, formula=PROD~IIND, family="gaussian")
summary(C7)
C8 <- glm(data=data.GLM02, formula=PROD~SHDI, family="gaussian")
summary(C8)

# MODEL SELECTION

# Univ models selection
UVAR.SLCT2 <- model.sel(C0, C1, C2, C5, C6, C7)

UVAR.SLCT2 

write.table(UVAR.SLCT2, "Data/UnivModelsCity.csv")

# # MULTIV MODELS
# CA <- glm(data=data.GLM02, formula=PROD~IIND*CVRG, family="gaussian")
# summary(CA)
# CB <- glm(data=data.GLM02, formula=PROD~IIND*DNST, family="gaussian")
# summary(CB)
# CC <- glm(data=data.GLM02, formula=PROD~IIND*PRXM, family="gaussian")
# summary(CC)
# CD <- glm(data=data.GLM02, formula=PROD~IIND*CNNC, family="gaussian")
# summary(CD)
# CE <- glm(data=data.GLM02, formula=PROD~IIND*SHEI, family="gaussian")
# summary(CE)
# CF <- glm(data=data.GLM02, formula=PROD~CVRG*SHEI, family="gaussian")
# summary(CF)
# 
# # Multiv models selection
# MVAR.SLCT2 <- model.sel(C0, CB, CC, CD, CE)
# 
# MVAR.SLCT2
# 
# write.table(MVAR.SLCT2, "MulvModelsCity.csv")

############################
### T-test national mean ###
############################
t.test(data.prop$ProdArea_kg, mu=1932)
t.test(data.muni$ProdArea_kg, mu=1932)

mean(data.prop$ProdArea_kg)
median(data.prop$ProdArea_kg)
sd(data.prop$ProdArea_kg)

mean(data.muni$ProdArea_kg)
median(data.muni$ProdArea_kg)
sd(data.muni$ProdArea_kg)
