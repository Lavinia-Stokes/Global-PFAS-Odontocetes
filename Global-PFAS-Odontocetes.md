# Global-PFAS-Odontocetes

##### Total Body-length Index #####

#### (as a proxy to estimate life stage of odontocetes) ####


# Load necessary libraries
library(dplyr)
library(scales)  # for rescale function

# Load individual dolphin dataset
sorted_data <- read.csv("sorted_data.csv")

# Clean 'Sex' column
sorted_data <- sorted_data %>%
  mutate(Sex = case_when(
    Sex %in% c("TBC", "", "0") ~ NA_character_,
    TRUE ~ Sex
  ))

# Convert 'Size' column to numeric and replace 0 values with NA
sorted_data <- sorted_data %>%
  mutate(Size = as.numeric(as.character(Size)),
         Size = ifelse(Size == 0, NA, Size))

# Filter relevant data
subset_data <- sorted_data %>%
  filter(!is.na(Size) & Size != "" & Species != "")

# Load species size reference dataset
ldata <- read.csv("ldata_revised.csv")

# Clean column names
colnames(ldata) <- gsub(" ", "_", colnames(ldata))  # Remove spaces
colnames(ldata) <- trimws(colnames(ldata))          # Trim hidden characters
colnames(ldata) <- make.names(colnames(ldata))      # Make names safe for R

# Define expected columns
expected_cols <- c("Species", "Min", "F_Max", "M_Max")

# Check for missing columns
available_cols <- intersect(expected_cols, colnames(ldata))
missing_cols <- setdiff(expected_cols, available_cols)

if (length(missing_cols) > 0) {
  print("ERROR: The following expected columns are missing in ldata:")
  print(missing_cols)
  stop("Check the printed column names and rename them if needed!")
}

# Select only required columns
ldata <- ldata[, available_cols]

# Ensure size columns are numeric
ldata <- ldata %>%
  mutate(Min = as.numeric(Min),
         F_Max = as.numeric(F_Max),
         M_Max = as.numeric(M_Max))

# Merge datasets on species
merged_data <- subset_data %>%
  left_join(ldata, by = "Species")

# Initialize the INDEX column
merged_data$INDEX <- NA

# Compute index for females (manually scaling)
merged_data <- merged_data %>%
  mutate(INDEX = ifelse(Sex == "F" & !is.na(F_Max) & !is.na(Min),
                        (Size - Min) / (F_Max - Min),
                        INDEX))

# Compute index for males (manually scaling)
merged_data <- merged_data %>%
  mutate(INDEX = ifelse(Sex == "M" & !is.na(M_Max) & !is.na(Min),
                        (Size - Min) / (M_Max - Min),
                        INDEX))

# Handle unknown sex by averaging male and female indices
unknown_sex_samples <- c("D_99", "D_46", "D_83", "D_96", "D_119", "D_129")

merged_data <- merged_data %>%
  group_by(Species) %>%
  mutate(
    avg_index = mean(c(INDEX[Sex == "F"], INDEX[Sex == "M"]), na.rm = TRUE),
    INDEX = ifelse(Sample %in% unknown_sex_samples, avg_index, INDEX)
  ) %>%
  ungroup()

# Cap index values between 0 and 1 to ensure valid scaling
merged_data$INDEX <- pmax(0, pmin(merged_data$INDEX, 1))

# Remove unwanted samples
merged_data <- merged_data %>%
  filter(!Sample %in% c("D_772", "D_880") & !is.na(INDEX))

# Write the final dataset to CSV
write.csv(merged_data, file = "PFAS_index_final.csv", row.names = FALSE)

# View relevant columns
index_view <- dplyr::select(merged_data, Paper, Sex, Species, Size, Sample, Common.name, INDEX)
print(index_view)


# Ensure size columns are numeric (avoiding character issues)
merged_data <- merged_data %>%
  mutate(Size = as.numeric(Size),
         Min = as.numeric(Min),
         F_Max = as.numeric(F_Max),
         M_Max = as.numeric(M_Max))
# Check if any Min or Max values are missing
print("Missing values in Min column:")
print(sum(is.na(merged_data$Min)))

print("Missing values in F_Max column:")
print(sum(is.na(merged_data$F_Max)))

print("Missing values in M_Max column:")
print(sum(is.na(merged_data$M_Max)))
# Compute index for females
merged_data <- merged_data %>%
  mutate(INDEX = ifelse(Sex == "F" & !is.na(F_Max) & !is.na(Min),
                        (Size - Min) / (F_Max - Min),
                        INDEX))

# Compute index for males
merged_data <- merged_data %>%
  mutate(INDEX = ifelse(Sex == "M" & !is.na(M_Max) & !is.na(Min),
                        (Size - Min) / (M_Max - Min),
                        INDEX))

# Check if any values are out of range
print("Index values less than 0:")
print(sum(merged_data$INDEX < 0, na.rm = TRUE))

print("Index values greater than 1:")
print(sum(merged_data$INDEX > 1, na.rm = TRUE))

# Cap index values between 0 and 1
merged_data$INDEX <- pmax(0, pmin(merged_data$INDEX, 1))



# glmm-pfas-odontocetes

####### sumPFAS analysis #########

#Set working directory
setwd("/Users/laviniastokes/Documents/A HONS/Global Analysis/Manuscript")
# Remove all objects from the environment to ensure a clean workspace
rm(list = ls())

# Load required libraries
library(lme4)
library(glmmTMB)
library(performance)
library(ggeffects)
library(ggplot2)
library(AICcmodavg)
library(dplyr)
#library(model)
library(flexplot)

# Load cleaned data
dat <- read.csv("cleandata_JUNE.csv", header = TRUE, sep = ",")

# Convert categorical variables to factors
dat$Group <- as.factor(dat$Group)
dat$Sex <- as.factor(dat$Sex)
dat$Ocean_loc <- as.factor(dat$Ocean_loc)

# Standardize continuous predictors
dat$INDEX_scaled <- scale(dat$INDEX)
dat$Year_scaled <- scale(dat$Year)

# Filter out rows with missing values in any relevant variables
dat_clean <- dat %>% 
  filter(!is.na(totalPFAS), !is.na(Sex), !is.na(INDEX_scaled), !is.na(Year_scaled), !is.na(Ocean_loc))

colnames(dat_clean)<-c("ref","species","latitude","longitude","sampleID","group","sex","index","ocean","country","year","PFCAs","PFSAs","totalPFAS","index_scaled","year_scaled")

# If your data are sorted by time, use acf (Autocorrelation Function)
dat_sorted <- dat_clean %>% arrange(year)
head(dat_sorted[c("year", "totalPFAS")])

# Now compute the autocorrelation function (ACF)
out_ACF<-acf(dat_sorted$totalPFAS, main = "ACF of totalPFAS (sorted by year)")
# Given this result, the PFAS values across years do not show strong autocorrelation
# best option Use year_scaled as a fixed effect in the model to:
#.    + Account for long-term trend in PFAS (if any),
#.    + Prevent potential confounding due to year-to-year variation.
# the ACF shows me that PFAS values across years are not serially correlated, meaning they don’t "depend" on the value from the previous year.


##==================================================================================
## STEP #1: COMPUTING A BASELINE MODEL (RANDOM EFFECT ANOVA) = GROUP & LOCATION
##==================================================================================
#Creating this baseline is useful because it allows us to ask:
#How much of the variation in PFAS concentrations is explained just
#by species group or location, without any other covariates?
#You can compare more complex models (with fixed effects like Year, Sex, Age Index)
#against this baseline using metrics like AIC, BIC, or likelihood ratio tests.
#to help evaluate whether adding predictors improves the model fit.

### GROUPS ##########
hist(dat_clean$totalPFAS, breaks = 50, main = "Distribution of totalPFAS", xlab = "totalPFAS")
mean(dat_clean$totalPFAS == 0)  # proportion  
boxplot(totalPFAS ~ group, data = dat_clean, main = "Boxplot of totalPFAS by Group", xlab = "Group", ylab = "totalPFAS")

#flip the axes so we can see the spread of PFAS across all groups
ggplot(dat_clean, aes(x = group, y = totalPFAS)) +
  geom_boxplot() +
  coord_flip() +  # flips the x and y axes
  theme_minimal() +
  labs(x = "Group", y = "total PFAS", title = "Boxplot of total PFAS by Group")

#The histograms and boxplots that are being plotted:
#Show the distribution and skew of total PFAS 
#(log-transform might be needed if it's right-skewed).
#Show variation between groups before modeling.

#which can help: Decide on a transformation (e.g., log PFAS).
#Spot outliers or zero-inflated values.
#See if Group or Location have visibly different patterns.


# Fit models using the same dataset
#### Use Tweedie distribution instead of Gamma or Gaussian (log-link) as all Zeroes are included now - 
#### Tweedie allows exact zeroes while also modelling positive concentrations
base <- glmmTMB(totalPFAS ~ 1 + (1 | group), 
                family = tweedie(link = "log"), data = dat_clean)

# Check the model's fit and assumptions
check_model(base)

###### MODEL CHECK RESULTS GROUPS - my interpretations:
#Model captures skewed distribution of totalPFAS, Tweedie dist is good to use 
#Homogeneity: quite a bit of heteroscedasticity
#Uniformity of Residuals: Mostly follows linear line - seems to under-predict very low values (zeroes)
#Normality of Random Effects: Points are quite close to line, random effects dist. looks okay?

### LOCATION (OCEAN)######
# Check Location counts:
table(dat_clean$ocean) #Looks pretty good - Smallest group is South Atlantic (12), largest group is Pacific (227)

# Relevel so Oceania becomes the reference - 
# Reasons for this decision: Highest PFAS levels known within this group (Aus)
dat_clean$ocean <- relevel(dat_clean$ocean, ref = "Oceania")
base2 <- glmmTMB(totalPFAS ~ 1 + (1 | ocean), 
                 family = tweedie(link = "log"), data = dat_clean)

#	Ocean-level variation in PFAS exists (Std.Dev = 1.167), justifying its inclusion
exp(5.08)
#Average log PFAS ≈ 5.08 → ~161 ng/g on natural scale

##==================================================================================
## STEP #2: INTRACLASS CORRELATION (ICC) AND DESIGN EFFECT
##==================================================================================
## GROUP ########
# Calculate ICC to quantify the proportion of variance explained by random effects
#To answer the question:How much of the total variation in PFAS is explained by differences
#between groups (e.g., species or location), rather than within them?”.

var_comp <- VarCorr(base)$cond# Extract variance components

var_group <- attr(var_comp$group, "stddev")^2# Between-group variance (random intercept)

var_resid <- sigma(base)^2# Residual variance (for Tweedie, approximated from model output): 
#within-group (residual) variance — approximates differences between individuals within the same group.

# ICC = between-group variance / (between-group + residual)
ICC <- var_group / (var_group + var_resid)
print(paste("ICC =", round(ICC, 3)))

#About 6.9% of the total variation in PFAS is explained by differences between species groups.

# Design Effect tells you how much the variance is inflated due to the grouping
# Average group size
avg_group_size <- mean(table(dat_clean$group))
DEFF <- 1 + (avg_group_size - 1) * ICC # Design Effect
print(paste("Design Effect =", round(DEFF, 3)))

# FS: 
# The intraclass correlation coefficient (ICC = 0.069) indicates that approximately 6.9% of the total variance in total PFAS concentrations 
# is attributable to differences between groups. While this ICC is moderate, the Design Effect (DEFF = 3.269) suggests that the grouping structure 
# (e.g., group-level clustering) has a substantial impact on standard errors, effectively inflating them more than threefold compared to a simple random sample. 
# This reinforces the need to account for group as a random effect in your GLMM to avoid biased inference and underestimated uncertainty.

## LOCATION (OCEAN) ########
# Check the model's fit and assumptions
check_model(base2)


# Calculate ICC to quantify the proportion of variance explained by random effects
var_comp2 <- VarCorr(base2)$cond# Extract variance components
var_group2 <- attr(var_comp2$ocean, "stddev")^2# Between-group variance (random intercept)
var_resid2 <- sigma(base2)^2# Residual variance (for Tweedie, approximated from model output)

# ICC = between-group variance / (between-group + residual)
ICC2 <- var_group2 / (var_group2 + var_resid2)

print(paste("ICC =", round(ICC2, 3)))

##About 3.5% of the total variation in PFAS is explained by differences between oceans. 

# Design Effect tells you how much the variance is inflated due to the grouping
# Average group size
avg_group_size2 <- mean(table(dat_clean$ocean))
DEFF2 <- 1 + (avg_group_size2 - 1) * ICC2 # Design Effect
print(paste("Design Effect =", round(DEFF2, 3)))

##FS = 
# ICC is lower than the species group ICC (0.069), meaning less between-ocean variance relative to total variance.
# On its own, this would suggest ocean contributes only modestly to explaining total PFAS variation.
# Despite the lower ICC, the large group sizes (especially Pacific: n = 227) dramatically inflate the design effect (5.143).
# This means that even small ICCs can cause major precision loss if group sizes are unbalanced or large.
# Including ocean as a random effect is statistically necessary, particularly to correct standard errors and avoid type I errors.
# It also helps control for spatial structure, making your fixed-effect estimates (e.g., for sex, index, year) more robust.

##==================================================================================
## STEP #3: MODEL FITS
##==================================================================================
## FIT1 = OCEAN AS RANDOM EFFECT
fit1 <- glmmTMB(totalPFAS ~ (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(base, fit1)# Compare full model with baseline model

Cand.set <- list(base, fit1)  # List of candidate models
Modnames <- c("Base Model", "Fit1")  # Model names
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out, model.high = "Fit1", model.low = "Base Model")

##FS = 
# Fit1 (with both ocean and group as random effects) is over 940,000 times more supported than the base model with only group.
# The ΔAICc = 27.5 and AICc weight = 1.0 confirm this is a substantial and statistically compelling improvement.
# This is driven by: (i) The large design effect from ocean (DEFF > 5), and (ii) the extra variance captured by adding ocean as a random effect

check_model(fit1)   # Visualize model diagnostics and fixed effect estimates
summary(fit1)

## FIT2 
## = fixes effects = SEX 
## = random effects = OCEAN and Group
fit2 <- glmmTMB(totalPFAS ~ sex + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit1, fit2)# Compare full model with baseline model

Cand.set <- list(base, fit1, fit2)  # List of candidate models
Modnames <- c("Base Model", "Fit1", "Fit2")  # Model names
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out, model.high = "Fit2", model.low = "Fit1")

#FS = Fit2 (which adds sex as a fixed effect) is 3.7× more supported than Fit1.
#     The effect of sexM is indicates that males have higher PFAS concentrations. 
#     Estimate = 0.224 (log scale) → On average, males have exp(0.224) ≈ 1.25× higher total PFAS than females.

check_model(fit2)   # Visualize model diagnostics and fixed effect estimates
estimates(fit2) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit2)
#FS = Minor heteroscedasticity at low values, which is expected and acceptable under the Tweedie distribution.
#.    Model fit is visually appropriate, with residuals mostly well-behaved.


## FIT3 
## = fixes effects = SEX & Index
## = random effects = OCEAN and Group
fit3 <- glmmTMB(totalPFAS ~ sex + index_scaled + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit2, fit3)# Compare full model with baseline model

Cand.set <- list(base, fit1, fit2, fit3)  # List of candidate models
Modnames <- c("Base Model", "Fit1", "Fit2","Fit3")  # Model names
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out, model.high = "Fit3", model.low = "Fit2")
#FS =  Fit3 significantly improves model fit: it is 1,279× more supported than Fit2, with ΔAICc = 14.3.
#     Fit3 is over a thousand times more likely to be the best model compared to Fit2.
#     Both sex and index_scaled are statistically significant predictors.
#     The effect of index_scaled is negative, indicating higher ecological index values are associated with lower PFAS.

check_model(fit3)   # Visualize model diagnostics and fixed effect estimates
estimates(fit3) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit3)

#FS = Sex and ecological index are both meaningful predictors of PFAS.
#     Higher index (possibly indicating better condition, less stress, or ecological richness) correlates with lower PFAS levels, 
#     suggesting possible behavioral or habitat-mediated exposure differences.
#     sexM = 0.235, p = 0.022 → Male individuals have ~26% higher PFAS (exp(0.235) ≈ 1.26×) than females.
#     index_scaled = -0.171, p < 0.001 → For each 1 SD increase in index, PFAS decreases by ~16% (exp(-0.171) ≈ 0.84×).


##### FIT4 
## = fixes effects = SEX & Index & Year
## = random effects = OCEAN and Group
fit4 <- glmmTMB(totalPFAS ~ sex + index_scaled + year_scaled + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit3, fit4)# Compare full model with baseline model



Cand.set <- list(base, fit1, fit2, fit3, fit4)  # List of candidate models
Modnames <- c("Base Model", "Fit1", "Fit2","Fit3", 'Fit4')  # Model names
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out, model.high = "Fit4", model.low = "Fit3")

check_model(fit4)   # Visualize model diagnostics and fixed effect estimates
estimates(fit4) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit4)

#Pseudo R²: Shows the marginal R² and conditional R² to show how much variance is explained by the random and fixed effects.
r2 <- r2_nakagawa(fit4)
print(r2)

r2_3 <- r2_nakagawa(fit3)
print(r2_3)

r2_5 <- r2_nakagawa(fit5)
print(r2_5)

r2_2 <- r2_nakagawa(fit2)
print(r2_2)

r2_1 <- r2_nakagawa(fit1)
print(r2_1)



## FS = The best-supported model (ΔAICc = 0; AICc weight = 0.73) 
#       males exhibited significantly higher PFAS concentrations than females (β = 0.251, p = 0.015), = 29% increase. 
#       Higher ecological index scores were associated with significantly lower PFAS levels (β = -0.162, p < 0.001), = ~15% reduction per standard deviation increase in the index. 
#       PFAS concentrations also showed a modest but significant increase over time (β = 0.129, p = 0.043). 
#       Random effects analysis confirmed substantial variance attributable to both group (SD = 1.47) and ocean (SD = 0.72), justifying their inclusion to correct for clustering effects. 
#       The model explained 76.1% of the total variance (conditional R²), with the fixed effects accounting for 1.4% (marginal R²), = much of the variation in PFAS burden is structured by group and spatial factors.

## FIT5 
## = fixes effects = SEX & Index & Year + Interaction (Sex * Index)
## = random effects = OCEAN and Group
fit5 <- glmmTMB(totalPFAS ~ sex + index_scaled + year_scaled + (sex*index_scaled) + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit4, fit5)# Compare full model with baseline model

Cand.set <- list(base, fit1, fit2, fit3, fit4, fit5)  # List of candidate models
Modnames <- c("Base Model", "Fit1", "Fit2","Fit3",'Fit4','Fit5')  # Model names
aicctable.out <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out)

length(unique(dat_clean$ref))



# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out, model.high = "Fit4", model.low = "Fit5")

# FS = Fit4 remains the best-supported model, with the lowest AICc.
#     Fit5 is 2.78× less supported than Fit4 (evidence ratio = 2.78).
#     The predicted difference range in Fit5 is small for most observations, suggesting limited added explanatory power from the interaction.Although adding the sex × index interaction (Fit5) slightly improves log-likelihood, it does not justify the extra complexity, as Fit4 has better AICc and higher weight.
#     Fit4 remains the most parsimonious model, and the interaction effect in Fit5 appears negligible for most of the distribution
#     conclusion = no need for the interaction


##==================================================================================
## STEP #4: FINAL VISUALIZATION AND CLEANUP
##==================================================================================
library(ggplot2)
ranef_fit4 <- ranef(fit4)

# Extract random effects
group_re <- ranef_fit4$cond$group
ocean_re <- ranef_fit4$cond$ocean

# Convert to data frames
group_df <- data.frame(group = rownames(group_re), effect = group_re[, 1])
ocean_df <- data.frame(ocean = rownames(ocean_re), effect = ocean_re[, 1])

# Plot: Group random effects
ggplot(group_df, aes(x = reorder(group, effect), y = effect)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Random Intercepts by Group", x = "Group", y = "Effect (log-scale deviation)")

# FS = strong variation in PFAS concentrations across cetacean taxa, after accounting for sex, ecological index, and year. 
#    Sousa, Tursiops, and Neophocaena = much higher-than-expected PFAS levels, with predicted concentrations approximately 9 to 18 times above the baseline. 
#    Lagenorhynchus, Pseudorca, and Orcinus = lower PFAS, with levels up to 80% below expected values. 
#    These differences suggest that species-specific ecological or physiological factors, such as foraging behavior or habitat use, play a significant role in PFAS exposure and accumulation. 
#    Emphasized on the importance of modeling group as a random effect to capture meaningful biological heterogeneity.

# Plot: Ocean random effects
ggplot(ocean_df, aes(x = reorder(ocean, effect), y = effect)) +
  geom_point(color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Random Intercepts by Ocean", x = "Ocean", y = "Effect (log-scale deviation)")

# Add a column to identify the random effect type
ocean_df$Type <- "Ocean"
group_df$Type <- "Genus"

# Rename the x-variable column to a common name, e.g. "Level"
ocean_df$Level <- ocean_df$ocean
group_df$Level <- group_df$group

# Combined Random Effects to plot together with two panels:

combined_df <- rbind(ocean_df[, c("Level", "effect", "Type")],
                     group_df[, c("Level", "effect", "Type")])


ggplot(combined_df, aes(x = reorder(Level, effect), y = effect)) +
  geom_point(color = "black", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ Type, scales = "free_y") +
  labs(
       x = NULL,
       y = "Effect (log-scale deviation)") +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10))
ggplot(combined_df, aes(x = reorder(Level, effect), y = effect, color = Type)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ Type, scales = "free_y") +
  scale_color_manual(values = c("Genus" = "black", "Ocean" = "#008ECC")) +
  labs(
    x = NULL,
    y = "Effect (log-scale deviation)"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10),
        legend.position = "none")  # hide legend if not needed

ggsave("random_effects_plot.pdf", width = 10, height = 6, dpi = 300)


# FS = spatial variation in PFAS concentrations after controlling for individual-level predictors. 
#      Pacific and Oceania show the highest positive deviations from the overall mean, with log-scale effects suggesting PFAS levels approximately 2.3× and 1.3× higher than average, respectively. 
#.     Similarly, the North Atlantic also shows elevated levels (~1.28×). 
#.     the Mediterranean has the strongest negative deviation, with PFAS concentrations estimated at only ~34% of the baseline, while the Arctic Ocean and South Atlantic exhibit values close to or slightly below the global average. 
#      These spatial patterns suggest that regional factors—such as industrial activity, oceanic currents, or pollutant transport pathways—may influence PFAS exposure = that is why this is importaant to consider ocean as a random effect.

## FINAL CONCLUSION
# Model = Best model among all candidates include = fixed effect = sex , index and year + random effect = ocean + group
#         76.1% of total variance explained (includes random effects) with Fixed effects alone explain only 1.4% of variance
#         Interpretation: Random effects dominate, indicating strong group- and region-specific variation in PFAS

# Fixed Effects
#         Males have ~29% higher PFAS concentrations than females (exp(0.251) ≈ 1.29)
#         Index (standardized): each 1 SD increase in index is associated with a ~15% decrease in PFAS (exp(-0.162) ≈ 0.85)
#         Year (standardized): ~14% increase in PFAS over time (exp(0.129) ≈ 1.14)

# Random Effects
#         Largest positive deviations: Sousa (~18× higher PFAS) > Tursiops (17x higher) > Neophocaena (9x higher)
#         Largest negative deviations: Lagenorhynchus (80%) > Pseudorca (73%) > Orcinus (74%) =  ~70–80% lower
#         Ocean Region: Pacific (~2.3× higher) > Oceania > North Atlantic  > South Atlantic > Arctic Ocean > Mediterranean (~66% lower)


# Create prediction data
PFAS_age_pred <- ggpredict(fit4, terms = "index_scaled")
PFAS_year_pred <- ggpredict(fit4, terms = "year_scaled")
PFAS_age_pred <- ggpredict(fit4, terms = "index_scaled", ci.lvl = 0.95)
PFAS_year_pred <- ggpredict(fit4, terms = "year_scaled", ci.lvl = 0.95)


# If your raw (unscaled) values are stored separately, e.g., in your dataset:
# You can back-transform manually if you know the original mean and SD:

mean_index <- mean(dat_clean$index, na.rm = TRUE)
sd_index <- sd(dat_clean$index, na.rm = TRUE)

mean_year <- mean(dat_clean$year, na.rm = TRUE)
sd_year <- sd(dat_clean$year, na.rm = TRUE)


# Back-transform x-values from scaled to original
PFAS_age_pred$x_unscaled <- PFAS_age_pred$x * sd_index + mean_index
PFAS_year_pred$x_unscaled <- PFAS_year_pred$x * sd_year + mean_year


# Age plot with log-scale PFAS
PFAS_ageplot <- ggplot(PFAS_age_pred, aes(x = x_unscaled, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#C5F6F7", alpha = 0.6) +
  geom_line(color = "black", size = 0.4) +
  labs(x = "Age Index", y = "∑PFAS (ng/g ww)") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.3),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    plot.tag = element_text(face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  labs(tag = "A")



# Year plot with log-scale PFAS
PFAS_yearplot <- ggplot(PFAS_year_pred, aes(x = x_unscaled, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#FECACA", alpha = 0.6) +
  geom_line(color = "black", size = 0.4) +
  labs(x = "Year", y = "∑PFAS (ng/g ww)") +
  theme_minimal() +
  theme(    panel.grid = element_blank(),
            axis.line = element_line(color = "black", linewidth =  0.3),
            axis.ticks = element_line(color = "black", linewidth =  0.3),
            plot.tag = element_text(face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10)
  ) +
  labs(tag = "B")

PFAS_ageplot + PFAS_yearplot  # Combine using patchwork or cowplot

# Combine plots (if you haven’t already)
combined_plot <- PFAS_ageplot / PFAS_yearplot

# Save to file
ggsave("PFAS_prediction_plots.png", combined_plot,
       width = 6, height = 7, dpi = 300, units = "in")

# Predict using unscaled terms
# (Assumes you used scaled variables in the model but want to plot unscaled)
PFAS_age_pred <- ggpredict(fit4, terms = "index_scaled [all]")  # get full range
PFAS_year_pred <- ggpredict(fit4, terms = "year_scaled [all]")


# Optional: Save to file
ggsave("log_PFAS_prediction_plots.png", combined_plot, width = 6, height = 7, dpi = 300)


#Results of total PFAS in a table
summary(fit4)
fit4_table <- tidy(fit4)
print(fit4_table)
write.csv(fit4_table, "fit4_summary_table.csv", row.names = FALSE)


###########Unique reference and species/ocean info for reference table##########

# Number of unique papers
length(unique(dat_clean$ref))
unique(dat_clean$ref)

# Number of unique genus groups
length(unique(dat_clean$group))
unique(dat_clean$group)

summary_by_genus <- dat_clean %>%
  group_by(group) %>%
  summarise(
    num_individuals = n(),                        # count rows (individuals)
    num_species = n_distinct(species)            # count unique species        # count unique ocean groups
  )
print(summary_by_genus)

# Number of unique species
length(unique(dat$Species))
unique(dat$Species)

#Number of unique countries
length(unique(dat_clean$country))
unique(dat_clean$country)

library(dplyr)

summary_by_paper <- dat_clean %>%
  group_by(ref) %>%
  summarise(
    num_individuals = n(),                        # count rows (individuals)
    num_species = n_distinct(species),            # count unique species
    num_ocean_groups = n_distinct(ocean)          # count unique ocean groups
  )

print(summary_by_paper)

species_ocean_by_paper <- dat_clean %>%
  group_by(ref) %>%
  summarise(
    species_list = paste(unique(species), collapse = ", "),
    ocean_groups = paste(unique(ocean), collapse = ", "),
    num_individuals = n()
  )

print(species_ocean_by_paper)

species_counts_per_paper <- dat_clean %>%
  group_by(ref, species) %>%
  summarise(num_individuals = n(), .groups = "drop")

print(species_counts_per_paper)

species_ocean_counts <- dat_clean %>%
  group_by(ref, ocean, species, country) %>%
  summarise(num_individuals = n(), .groups = "drop")

print(species_ocean_counts)





############PFCAs Vs PFOS##########


pfca_list <- c("PFCAs")  # adjust to your data
pfsa_list <- c("PFSAs")  # adjust to your data

class_summary <- dat_clean %>%
  mutate(Class = case_when(
    Compound %in% pfca_list ~ "PFCAs",
    Compound %in% pfsa_list ~ "PFSAs",
    TRUE ~ "Other"
  )) %>%
  filter(Class %in% c("PFCAs", "PFSAs")) %>%
  group_by(Class) %>%
  summarise(Total_Concentration = sum(Concentration, na.rm = TRUE))


library(ggplot2)

ggplot(dat_clean), aes(x = "PFCAs", "PFSAs", y = Total_Concentration, fill = Class)) +
  geom_col(width = 0.6) +
  labs(title = "Total PFAS Class Concentrations",
       x = "PFAS Class",
       y = "Total Concentration (ng/g)") +
  theme_minimal() +
  theme(legend.position = "none") 
  
#################################


######### Sum PFCAs analysis #############



#Set working directory
#setwd("/Users/laviniastokes/Documents/A HONS/Global Analysis/Manuscript")


# Load required libraries
library(lme4)
library(glmmTMB)
library(performance)
library(ggeffects)
library(ggplot2)
library(AICcmodavg)
library(dplyr)
#library(model)
library(flexplot)

# If your data are sorted by time, use acf (Autocorrelation Function)
dat_sorted <- dat_clean %>% arrange(year)
head(dat_sorted[c("year", "PFCAs")])

# Now compute the autocorrelation function (ACF)
out_ACF_PFCAs <-acf(dat_sorted$PFCAs, main = "ACF of PFCAs (sorted by year)")
# Given this result, your PFAS values across years do not show strong autocorrelation
# best option Use year_scaled as a fixed effect in your model to:
#.    + Account for long-term trend in PFAS (if any),
#.    + Prevent potential confounding due to year-to-year variation.


##==================================================================================
## STEP #1: COMPUTING A BASELINE MODEL (RANDOM EFFECT ANOVA) = GROUP & LOCATION
##==================================================================================
### GROUPS ##########
hist(dat_clean$PFCAs, breaks = 50, main = "Distribution of PFCAs", xlab = "PFCAs")
mean(dat_clean$PFCAs == 0)  # proportion  
boxplot(PFCAs ~ group, data = dat_clean, main = "Boxplot of PFCAs by Group", xlab = "Group", ylab = "PFCAs")

# Fit models using the same dataset
#### Use Tweedie distribution instead of Gamma or Gaussian (log-link) as all Zeroes are included now - 
#### Tweedie allows exact zeroes while also modelling positive concentrations
base_PFCAs <- glmmTMB(PFCAs ~ 1 + (1 | group), 
                family = tweedie(link = "log"), data = dat_clean)

# Check the model's fit and assumptions
check_model(base_PFCAs)

####### MODEL CHECK RESULTS GROUPS - my interpretations:
#Model captures skewed distribution of totalPFAS, Tweedie dist is good to use 
#Homogeneity: quite a bit of heteroscedasticity
#Uniformity of Residuals: Mostly follows linear line - seems to under-predict very low values (zeroes) and overpredict some values.
#Normality of Random Effects: Points are quite close to line, random effects dist. looks okay?

### LOCATION (OCEAN)##########
# Check Location counts:
table(dat_clean$ocean) #Looks pretty good - Smallest group is South Atlantic (12), largest group is Pacific (227)

# Relevel so Oceania becomes the reference - 
# Reasons for this decision: Highest PFAS levels known within this group (Aus)
dat_clean$ocean <- relevel(dat_clean$ocean, ref = "Oceania")
base2_PFCAs <- glmmTMB(PFCAs ~ 1 + (1 | ocean), 
                 family = tweedie(link = "log"), data = dat_clean)


##==================================================================================
## STEP #2: INTRACLASS CORRELATION (ICC) AND DESIGN EFFECT
##==================================================================================
## GROUP ########
# Calculate ICC to quantify the proportion of variance explained by random effects
var_comp_PFCAs <- VarCorr(base_PFCAs)$cond# Extract variance components
var_group_PFCAs <- attr(var_comp_PFCAs$group, "stddev")^2# Between-group variance (random intercept)
var_resid_PFCAs <- sigma(base_PFCAs)^2# Residual variance (for Tweedie, approximated from model output)

# ICC = between-group variance / (between-group + residual)
ICC_PFCAs <- var_group_PFCAs / (var_group_PFCAs + var_resid_PFCAs)
print(paste("ICC_PFCAs =", round(ICC_PFCAs, 3)))

# Design Effect tells you how much the variance is inflated due to the grouping
# Average group size
avg_group_size <- mean(table(dat_clean$group))
DEFF_PFCAs <- 1 + (avg_group_size - 1) * ICC_PFCAs # Design Effect
print(paste("Design Effect PFCAs =", round(DEFF_PFCAs, 3)))

# FS: 
# The intraclass correlation coefficient (ICC = 0.021) indicates that approximately 2.1% of the total variance in total PFAS concentrations 
# is attributable to differences between groups. While this ICC is moderate, the Design Effect (DEFF = 1.69) suggests that the grouping structure 
# (e.g., group-level clustering) has a substantial impact on standard errors, effectively inflating them around 1.3 times ( standard errors are ~√1.69 ≈ 1.3 times larger) 
#compared to a simple random sample. 
# This reinforces the need to account for group as a random effect in your GLMM to avoid biased inference and underestimated uncertainty.

## LOCATION (OCEAN) ########
# Check the model's fit and assumptions
check_model(base2_PFCAs)
# Calculate ICC to quantify the proportion of variance explained by random effects
var_comp2_PFCAs <- VarCorr(base2_PFCAs)$cond# Extract variance components
var_group2_PFCAs <- attr(var_comp2_PFCAs$ocean, "stddev")^2# Between-group variance (random intercept)
var_resid2_PFCAs <- sigma(base2_PFCAs)^2# Residual variance (for Tweedie, approximated from model output)

# ICC = between-group variance / (between-group + residual)
ICC2_PFCAs <- var_group2_PFCAs / (var_group2_PFCAs + var_resid2_PFCAs)
print(paste("ICC PFCAs =", round(ICC2_PFCAs, 3)))

# Design Effect tells you how much the variance is inflated due to the grouping
# Average group size
avg_group_size2_PFCAs <- mean(table(dat_clean$ocean))
DEFF2_PFCAs <- 1 + (avg_group_size2_PFCAs - 1) * ICC2_PFCAs # Design Effect
print(paste("Design Effect PFCAs =", round(DEFF2_PFCAs, 3)))

##FS = 
# ICC is lower than the group ICC (0.013), meaning less between-ocean variance relative to total PFCAs variance.
# On its own, this would suggest ocean contributes only modestly to explaining total PFCAs variation.
# Despite the lower ICC, the large group sizes (especially Pacific: n = 227) inflate the design effect (2.564).
# This means that even small ICCs can cause major precision loss if group sizes are unbalanced or large.
# Including ocean as a random effect is statistically necessary, particularly to correct standard errors and avoid type I errors.
# It also helps control for spatial structure, making your fixed-effect estimates (e.g., for sex, index, year) more robust.

##==================================================================================
## STEP #3: MODEL FITS
##==================================================================================
## FIT1 = OCEAN AS RANDOM EFFECT
fit1_PFCAs <- glmmTMB(PFCAs ~ (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(base_PFCAs, fit1_PFCAs)# Compare full model with baseline model

Cand.set_PFCAs <- list(base_PFCAs, fit1_PFCAs)  # List of candidate models
Modnames_PFCAs <- c("Base Model PFCAs", "Fit1 PFCAs")  # Model names
aicctable.out_PFCAs <- aictab(cand.set = Cand.set_PFCAs, modnames = Modnames_PFCAs)  # AICc comparison
print(aicctable.out_PFCAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFCAs, model.high = "Fit1 PFCAs", model.low = "Base Model PFCAs")

##FS = 
# Fit1 (with both ocean and group as random effects) is over 3.87 trillion times more supported than the base model with only group.
# The ΔAICc = 57.97 and AICc weight = 1.00 confirms this is a substantial and statistically compelling improvement.
# This is driven by: (i) The large design effect from ocean (DEFF > 2), and (ii) the extra variance captured by adding ocean as a random effect

check_model(fit1_PFCAs)   # Visualize model diagnostics and fixed effect estimates
summary(fit1_PFCAs)

## FIT2 
## = fixes effects = SEX 
## = random effects = OCEAN and Group
fit2_PFCAs <- glmmTMB(PFCAs ~ sex + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit1_PFCAs, fit2_PFCAs)# Compare full model with baseline model

Cand.set <- list(base_PFCAs, fit1_PFCAs, fit2_PFCAs)  # List of candidate models
Modnames <- c("Base Model PFCAs", "Fit1 PFCAs", "Fit2 PFCAs")  # Model names
aicctable.out_PFCAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFCAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFCAs, model.high = "Fit2 PFCAs", model.low = "Fit1 PFCAs")

#FS = Fit1 is 1.9× (1 / 0.52 ≈ 1.92) more supported than Fit2.
# Adding sex does not substantially improve model fit 

check_model(fit2_PFCAs)   # Visualize model diagnostics and fixed effect estimates
estimates(fit2_PFCAs) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit2_PFCAs)

#FS = Minor heteroscedasticity at low values, which is expected and acceptable under the Tweedie distribution.
#.    Model fit is visually appropriate, with residuals mostly well-behaved.


## FIT3 
## = fixed effects = SEX & Index
## = random effects = OCEAN and Group
fit3_PFCAs <- glmmTMB(PFCAs ~ index_scaled + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit2_PFCAs, fit3_PFCAs)# Compare full model with baseline model

Cand.set <- list(base_PFCAs, fit1_PFCAs, fit2_PFCAs, fit3_PFCAs)  # List of candidate models
Modnames <- c("Base Model", "Fit1 PFCAs", "Fit2 PFCAs","Fit3 PFCAs")  # Model names
aicctable.out_PFCAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFCAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFCAs, model.high = "Fit3 PFCAs", model.low = "Fit1 PFCAs")
#FS =  Fit3 significantly improves model fit: it is 81537724× more supported than Fit1, with ΔAICc = 36.43
#     Fit3 is way more likely to be the best model compared to Fit2.
#     index_scaled is a statistically significant predictor.
#     The effect of index_scaled is negative, indicating higher ecological index values are associated with lower PFAS.

check_model(fit3_PFCAs)   # Visualize model diagnostics and fixed effect estimates
estimates(fit3_PFCAs) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit3_PFCAs)

#FS = Ecological index is a meaningful predictor of PFAS.
#     Higher index (possibly indicating better condition, less stress, or ecological richness) correlates with lower PFAS levels, 
#     suggesting possible behavioral or habitat-mediated exposure differences.
#     index_scaled = -0.23420, p < 0 → For each 1 SD increase in index, PFAS decreases by ~20.9% (exp(-0.2342) ≈ 0.79×).

##### FIT4 
## = fixes effects = SEX & Index & Year
## = random effects = OCEAN and Group
fit4_PFCAs <- glmmTMB(PFCAs ~ index_scaled + year_scaled + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit3_PFCAs, fit4_PFCAs)# Compare full model with baseline model



Cand.set <- list(base_PFCAs, fit1_PFCAs, fit2_PFCAs, fit3_PFCAs, fit4_PFCAs)  # List of candidate models
Modnames <- c("Base Model PFCAs", "Fit1 PFCAs", "Fit2 PFCAs","Fit3 PFCAs", "Fit4 PFCAs")  # Model names
aicctable.out_PFCAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFCAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFCAs, model.high = "Fit4 PFCAs", model.low = "Fit3 PFCAs")

check_model(fit4_PFCAs)   # Visualize model diagnostics and fixed effect estimates
estimates(fit4_PFCAs) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit4_PFCAs)

#Pseudo R²: Shows the marginal R² and conditional R² to show how much variance is explained by the random and fixed effects.
r2_4_PFCAs <- r2_nakagawa(fit4_PFCAs)
print(r2_4_PFCAs)

#Random effects explain approx 74.6% of the total variation of PFCAs (Conditional = 0.746)
#Fixed effects explain approx 5.7% of the total variation of PFCAs (Marginal = 0.057)

r2_3_PFCAs <- r2_nakagawa(fit3_PFCAs)
print(r2_3_PFCAs)

r2_5_PFCAs <- r2_nakagawa(fit5_PFCAs)
print(r2_5_PFCAs)

r2_2_PFCAs <- r2_nakagawa(fit2_PFCAs)
print(r2_2_PFCAs)

r2_1_PFCAs <- r2_nakagawa(fit1_PFCAs)
print(r2_1_PFCAs)


## FS = The best-supported model (ΔAICc = 0; AICc weight = 1) 
#       Higher ecological index scores were associated with significantly lower PFAS levels (β = -0.18808, p < 0), = ~17.1% reduction per standard deviation increase in the index. 
#       PFAS concentrations also showed a significant increase over time (β = 0.37070, p < 0), = ~ 44.9% increaseof PFCAs per unit of time
#       Random effects analysis confirmed substantial variance attributable to both group (SD = 1.2044) and ocean (SD = 0.8457), justifying their inclusion to correct for clustering effects. 
#       The model explained 74.6% of the total variance (conditional R²), with the fixed effects accounting for 5.7% (marginal R²), = much of the variation in PFAS burden is structured by group and spatial factors.

##==================================================================================
## STEP #4: FINAL VISUALIZATION AND CLEANUP
##==================================================================================
library(ggplot2)
ranef_fit4_PFCAs <- ranef(fit4_PFCAs)

# Extract random effects
group_re_PFCAs <- ranef_fit4_PFCAs$cond$group

# Sousa 6x higher, Neophocaena 5x, Tursiops 3x
# Pontoporia 71% lower, Lagenorhynchys 82%, Cephalorhynchus 83%

ocean_re_PFCAs <- ranef_fit4_PFCAs$cond$ocean
exp(1.1062113)

# Higher: Pacific 3x, North Atlantic 1.5x, Mediterranean 1.4x
#Lower: South Atlantic 50% below mean, Oceania 67% below mean

# Convert to data frames
group_df_PFCAs <- data.frame(group = rownames(group_re_PFCAs), effect = group_re_PFCAs[, 1])
ocean_df_PFCAs <- data.frame(ocean = rownames(ocean_re_PFCAs), effect = ocean_re_PFCAs[, 1])

# Plot: Group random effects
ggplot(group_df_PFCAs, aes(x = reorder(group, effect), y = effect)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Random Intercepts by Group", x = "Group", y = "Effect (log-scale deviation)")

# FS = strong variation in PFAS concentrations across cetacean taxa, after accounting for sex, ecological index, and year. 
#    Sousa, Tursiops, and Neophocaena = much higher-than-expected PFAS levels, with predicted concentrations approximately 9 to 18 times above the baseline. 
#    Lagenorhynchus, Pseudorca, and Orcinus = lower PFAS, with levels up to 80% below expected values. 
#    These differences suggest that species-specific ecological or physiological factors, such as foraging behavior or habitat use, play a significant role in PFAS exposure and accumulation. 
#    Emphasized on the importance of modeling group as a random effect to capture meaningful biological heterogeneity.

# Plot: Ocean random effects
ggplot(ocean_df_PFCAs, aes(x = reorder(ocean, effect), y = effect)) +
  geom_point(color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Random Intercepts by Ocean", x = "Ocean", y = "Effect (log-scale deviation)")

# Add a column to identify the random effect type
ocean_df_PFCAs$Type <- "Ocean"
group_df_PFCAs$Type <- "Genus"

# Rename the x-variable column to a common name, e.g. "Level"
ocean_df_PFCAs$Level <- ocean_df_PFCAs$ocean
group_df_PFCAs$Level <- group_df_PFCAs$group

# Combined Random Effects to plot together with two panels:

combined_df_PFCAs <- rbind(ocean_df_PFCAs[, c("Level", "effect", "Type")],
                     group_df_PFCAs[, c("Level", "effect", "Type")])


ggplot(combined_df_PFCAs, aes(x = reorder(Level, effect), y = effect)) +
  geom_point(color = "black", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ Type, scales = "free_y") +
  labs(
    x = NULL,
    y = "Effect (log-scale deviation)") +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10))

ggsave("random_effects_plot_PFCAs.pdf", width = 10, height = 6, dpi = 300)


#################################


######### PFOS Analysis #############


#Set working directory
#setwd("/Users/laviniastokes/Documents/A HONS/Global Analysis/Manuscript")


# Load required libraries
library(lme4)
library(glmmTMB)
library(performance)
library(ggeffects)
library(ggplot2)
library(AICcmodavg)
library(dplyr)
#library(model)
library(flexplot)

# If your data are sorted by time, use acf (Autocorrelation Function)
dat_sorted <- dat_clean %>% arrange(year)
head(dat_sorted[c("year", "PFSAs")])

# Now compute the autocorrelation function (ACF)
out_ACF_PFSAs <-acf(dat_sorted$PFSAs, main = "ACF of PFSAs (sorted by year)")
# Given this result, your PFAS values across years do not show strong autocorrelation
# best option Use year_scaled as a fixed effect in your model to:
#.    + Account for long-term trend in PFAS (if any),
#.    + Prevent potential confounding due to year-to-year variation.


##==================================================================================
## STEP #1: COMPUTING A BASELINE MODEL (RANDOM EFFECT ANOVA) = GROUP & LOCATION
##==================================================================================
### GROUPS ##########
hist(dat_clean$PFSAs, breaks = 50, main = "Distribution of PFSAs", xlab = "PFSAs")
mean(dat_clean$PFSAs == 0)  # proportion  
boxplot(PFSAs ~ group, data = dat_clean, main = "Boxplot of PFSAs by Group", xlab = "Group", ylab = "PFSAs")

# Fit models using the same dataset
#### Use Tweedie distribution instead of Gamma or Gaussian (log-link) as all Zeroes are included now - 
#### Tweedie allows exact zeroes while also modelling positive concentrations
base_PFSAs <- glmmTMB(PFSAs ~ 1 + (1 | group), 
                      family = tweedie(link = "log"), data = dat_clean)

# Check the model's fit and assumptions
check_model(base_PFSAs)

####### MODEL CHECK RESULTS GROUPS - my interpretations:
#Model captures skewed distribution of totalPFAS, Tweedie dist is good to use 
#Homogeneity: quite a bit of heteroscedasticity
#Uniformity of Residuals: Mostly follows linear line - seems to under-predict very low values (zeroes) and overpredict some values.
#Normality of Random Effects: Points are quite close to line, random effects dist. looks okay?

### LOCATION (OCEAN)##########
# Check Location counts:
table(dat_clean$ocean) #Looks pretty good - Smallest group is South Atlantic (12), largest group is Pacific (227)

# Relevel so Oceania becomes the reference - 
# Reasons for this decision: Highest PFAS levels known within this group (Aus)
dat_clean$ocean <- relevel(dat_clean$ocean, ref = "Oceania")
base2_PFSAs <- glmmTMB(PFSAs ~ 1 + (1 | ocean), 
                       family = tweedie(link = "log"), data = dat_clean)


##==================================================================================
## STEP #2: INTRACLASS CORRELATION (ICC) AND DESIGN EFFECT
##==================================================================================
## GROUP ########
# Calculate ICC to quantify the proportion of variance explained by random effects
var_comp_PFSAs <- VarCorr(base_PFSAs)$cond# Extract variance components
var_group_PFSAs <- attr(var_comp_PFSAs$group, "stddev")^2# Between-group variance (random intercept)
var_resid_PFSAs <- sigma(base_PFSAs)^2# Residual variance (for Tweedie, approximated from model output)

# ICC = between-group variance / (between-group + residual)
ICC_PFSAs <- var_group_PFSAs / (var_group_PFSAs + var_resid_PFSAs)
print(paste("ICC_PFSAs =", round(ICC_PFSAs, 3)))

# Design Effect tells you how much the variance is inflated due to the grouping
# Average group size
avg_group_size <- mean(table(dat_clean$group))
DEFF_PFSAs <- 1 + (avg_group_size - 1) * ICC_PFSAs # Design Effect
print(paste("Design Effect PFSAs =", round(DEFF_PFSAs, 3)))

# FS: 
# The intraclass correlation coefficient (ICC = 0.095) indicates that approximately 9.5% of the total variance in total PFOS concentrations 
# is attributable to differences between groups. While this ICC is moderate, the Design Effect (DEFF = 4.127) suggests that the grouping structure 
# (e.g., group-level clustering) has a substantial impact on standard errors, effectively inflating them around 2.03 times ( standard errors are ~√4.127 ≈ 2.03 times larger) 
#compared to a simple random sample. 
# This reinforces the need to account for group as a random effect in your GLMM to avoid biased inference and underestimated uncertainty.

## LOCATION (OCEAN) ########
# Check the model's fit and assumptions
check_model(base2_PFSAs)
# Calculate ICC to quantify the proportion of variance explained by random effects
var_comp2_PFSAs <- VarCorr(base2_PFSAs)$cond# Extract variance components
var_group2_PFSAs <- attr(var_comp2_PFSAs$ocean, "stddev")^2# Between-group variance (random intercept)
var_resid2_PFSAs <- sigma(base2_PFSAs)^2# Residual variance (for Tweedie, approximated from model output)

# ICC = between-group variance / (between-group + residual)
ICC2_PFSAs <- var_group2_PFSAs / (var_group2_PFSAs + var_resid2_PFSAs)
print(paste("ICC PFSAs =", round(ICC2_PFSAs, 3)))

# Design Effect tells you how much the variance is inflated due to the grouping
# Average group size
avg_group_size2_PFSAs <- mean(table(dat_clean$ocean))
DEFF2_PFSAs <- 1 + (avg_group_size2_PFSAs - 1) * ICC2_PFSAs # Design Effect
print(paste("Design Effect PFSAs =", round(DEFF2_PFSAs, 3)))

##FS = 
# ICC is lower than the group ICC (0.046), meaning less between-ocean variance relative to PFOS variance.
# On its own, this would suggest ocean contributes only modestly to explaining PFOS variation.
# Despite the lower ICC, the large group sizes (especially Pacific: n = 227) inflate the design effect (6.474).
# This means that even small ICCs can cause major precision loss if group sizes are unbalanced or large.
# Including ocean as a random effect is statistically necessary, particularly to correct standard errors and avoid type I errors.
# It also helps control for spatial structure, making your fixed-effect estimates (e.g., for sex, index, year) more robust.

##==================================================================================
## STEP #3: MODEL FITS
##==================================================================================
## FIT1 = OCEAN AS RANDOM EFFECT
fit1_PFSAs <- glmmTMB(PFSAs ~ (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(base_PFSAs, fit1_PFSAs)# Compare full model with baseline model

Cand.set_PFSAs <- list(base_PFSAs, fit1_PFSAs)  # List of candidate models
Modnames_PFSAs <- c("Base Model PFSAs", "Fit1 PFSAs")  # Model names
aicctable.out_PFSAs <- aictab(cand.set = Cand.set_PFSAs, modnames = Modnames_PFSAs)  # AICc comparison
print(aicctable.out_PFSAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFSAs, model.high = "Fit1 PFSAs", model.low = "Base Model PFSAs")

##FS = 
# Fit1 (with both ocean and group as random effects) is over 10790036995  times more supported than the base model with only group.
# The ΔAICc = 46.2 and AICc weight = 1.00 confirms this is a substantial and statistically compelling improvement.
# This is driven by: (i) The large design effect from ocean (DEFF > 2), and (ii) the extra variance captured by adding ocean as a random effect

check_model(fit1_PFSAs)   # Visualize model diagnostics and fixed effect estimates
summary(fit1_PFSAs)

## FIT2 
## = fixes effects = SEX 
## = random effects = OCEAN and Group
fit2_PFSAs <- glmmTMB(PFSAs ~ sex + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit1_PFSAs, fit2_PFSAs)# Compare full model with baseline model

Cand.set <- list(base_PFSAs, fit1_PFSAs, fit2_PFSAs)  # List of candidate models
Modnames <- c("Base Model PFSAs", "Fit1 PFSAs", "Fit2 PFSAs")  # Model names
aicctable.out_PFSAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFSAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFSAs, model.high = "Fit2 PFSAs", model.low = "Fit1 PFSAs")

#FS = RATIOS < 2 SUGGESTS THE MODELS HAVE SIMILAR SUPPORT FROM THE DATA.

check_model(fit2_PFSAs)   # Visualize model diagnostics and fixed effect estimates
estimates(fit2_PFSAs) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit2_PFSAs)

#FS = Minor heteroscedasticity at low values, which is expected and acceptable under the Tweedie distribution.
#.    Model fit is visually appropriate, with residuals mostly well-behaved.


## FIT3 
## = fixed effects = SEX & Index
## = random effects = OCEAN and Group
fit3_PFSAs <- glmmTMB(PFSAs ~ sex + index_scaled + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit2_PFSAs, fit3_PFSAs)# Compare full model with baseline model

Cand.set <- list(base_PFSAs, fit1_PFSAs, fit2_PFSAs, fit3_PFSAs)  # List of candidate models
Modnames <- c("Base Model", "Fit1 PFSAs", "Fit2 PFSAs","Fit3 PFSAs")  # Model names
aicctable.out_PFSAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFSAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFSAs, model.high = "Fit3 PFSAs", model.low = "Fit2 PFSAs")
#FS =  Fit3 is ~64× more supported than Fit2, with ΔAICc = 8.32
#     Fit3 is more likely to be the best model compared to Fit2.
#     index_scaled is a statistically significant predictor.
#     The effect of index_scaled is negative, indicating higher ecological index values are associated with lower PFAS.

check_model(fit3_PFSAs)   # Visualize model diagnostics and fixed effect estimates
estimates(fit3_PFSAs) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit3_PFSAs)

#FS = Ecological index is a meaningful predictor of PFOS.
#     Higher index (possibly indicating better condition, less stress, or ecological richness) correlates with lower PFOS levels, 
#     suggesting possible behavioral or habitat-mediated exposure differences.
#     Sex = 0.21511, p <0.05 -> Males have ~ 24% higher PFOS levels than females
#     index_scaled = -0.14295, p < 0.001 → For each 1 SD increase in index, PFOS decreases by ~13.3% (exp(-0.14295) ≈ ).

##### FIT4 
## = fixes effects = SEX & Index & Year
## = random effects = OCEAN and Group
fit4_PFSAs <- glmmTMB(PFSAs ~ sex + index_scaled + year_scaled + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit3_PFSAs, fit4_PFSAs)# Compare full model with baseline model



Cand.set <- list(base_PFSAs, fit1_PFSAs, fit2_PFSAs, fit3_PFSAs, fit4_PFSAs)  # List of candidate models
Modnames <- c("Base Model PFSAs", "Fit1 PFSAs", "Fit2 PFSAs","Fit3 PFSAs", "Fit4 PFSAs")  # Model names
aicctable.out_PFSAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFSAs)

# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFSAs, model.high = "Fit4 PFSAs", model.low = "Fit3 PFSAs")

check_model(fit4_PFSAs)   # Visualize model diagnostics and fixed effect estimates
estimates(fit4_PFSAs) #Looks pretty good? Some scattering in the Homogeneity of variance for lower values
summary(fit4_PFSAs)

## FIT5 
## = fixes effects = SEX & Index & Year + Interaction (Sex * Index)
## = random effects = OCEAN and Group
fit5_PFSAs <- glmmTMB(PFSAs ~ sex + index_scaled + (sex*index_scaled) + (1 | ocean) + (1 | group),family = tweedie(link = "log"),data = dat_clean)
model.comparison(fit3_PFSAs, fit5_PFSAs)# Compare full model with baseline model

Cand.set <- list(base_PFSAs, fit1_PFSAs, fit2_PFSAs, fit3_PFSAs, fit4_PFSAs, fit5_PFSAs)  # List of candidate models
Modnames <- c("Base Model PFSAs", "Fit1 PFSAs", "Fit2 PFSAs", "Fit3 PFSAs", "Fit4 PFSAs", "Fit5 PFSAs") # Model names
aicctable.out_PFSAs <- aictab(cand.set = Cand.set, modnames = Modnames)  # AICc comparison
print(aicctable.out_PFSAs)


# The evidence ratio can be interpreted as the number of times a given model is more parsimonious than a lower-ranked model.
evidence(aic.table = aicctable.out_PFSAs, model.high = "Fit4 PFSAs", model.low = "Fit5 PFSAs")

# FS = Fit 3 remains the best-supported model, with the lowest AICc.
#     Fit5 is 1.52× less supported than Fit4 (evidence ratio = 1.52).
#     The predicted difference range in Fit5 is small for most observations, suggesting limited added explanatory power from the interaction.Although adding the sex × index interaction (Fit5) slightly improves log-likelihood, it does not justify the extra complexity, as Fit4 has better AICc and higher weight.
#     Fit4 remains the most parsimonious model, and the interaction effect in Fit5 appears negligible for most of the distribution
#     conclusion = no need for the interaction

#Pseudo R²: Shows the marginal R² and conditional R² to show how much variance is explained by the random and fixed effects.
r2_3_PFSAs <- r2_nakagawa(fit3_PFSAs) 
print(r2_3_PFSAs)


#Random effects explain approx 81.4% of the total variation of PFOS (Conditional = 0.814)
#Fixed effects explain approx 0.7% of the total variation of PFOS (Marginal = 0.007)

r2_4_PFSAs <- r2_nakagawa(fit4_PFSAs)
print(r2_4_PFSAs)

r2_5_PFSAs <- r2_nakagawa(fit5_PFSAs)
print(r2_5_PFSAs)

r2_2_PFSAs <- r2_nakagawa(fit2_PFSAs)
print(r2_2_PFSAs)

r2_1_PFSAs <- r2_nakagawa(fit1_PFSAs)
print(r2_1_PFSAs)

summary(fit3_PFSAs)
## FS   FIT3 = The best-supported model (ΔAICc = 0; AICc weight = 0.44) 
#FS = Ecological index is a meaningful predictor of PFOS.
#     Higher index (possibly indicating better condition, less stress, or ecological richness) correlates with lower PFOS levels, 
#     suggesting possible behavioral or habitat-mediated exposure differences.
#     Sex = 0.21511, p <0.05 -> Males have ~ 24% higher PFOS levels than females
#     index_scaled = -0.14295, p < 0.001 → For each 1 SD increase in index, PFOS decreases by ~13.3% (exp(-0.14295) ≈ ).

#       Higher ecological index scores were associated with significantly lower PFOS levels (β = -0.14295, p < 0), = ~13.3% reduction per standard deviation increase in the index. 
#       Random effects analysis confirmed substantial variance attributable to both group (SD = 1.7989) and ocean (SD = 0.7322), justifying their inclusion to correct for clustering effects. 
#       The model explained 81.4% of the total variance (conditional R²), with the fixed effects accounting for 0.7% (marginal R²), = majority of the variation in PFAS burden is structured by group and spatial factors.

##==================================================================================
## STEP #4: FINAL VISUALIZATION AND CLEANUP
##==================================================================================
library(ggplot2)
ranef_fit4_PFSAs <- ranef(fit4_PFSAs)

# Extract random effects
group_re_PFSAs <- ranef_fit4_PFSAs$cond$group

exp(3.8525012) #Sousa: 47x
exp(3.0954700) #Neophocaena: 22x
exp(2.9952468) #Tursiops: 20x

exp(-1.8660471) #Pseudorca: 84% lower
exp(-1.6670562) #Beaked whale: 81% lower
exp(-1.5853125) #Lagenorhynchus: 79% lower

ocean_re_PFSAs <- ranef_fit4_PFSAs$cond$ocean

exp(0.8268470) #Oceania


exp(0.6994707) #Arctic


#Higher: Oceania 2.3x, Pacific 1.4x, North Atlantic 1.2x, South Atlantic slightly above mean
#Lower: Mediterranean 60%, Arctic 50%



# Convert to data frames
group_df_PFSAs <- data.frame(group = rownames(group_re_PFSAs), effect = group_re_PFSAs[, 1])
ocean_df_PFSAs <- data.frame(ocean = rownames(ocean_re_PFSAs), effect = ocean_re_PFSAs[, 1])

# Plot: Group random effects
ggplot(group_df_PFSAs, aes(x = reorder(group, effect), y = effect)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Random Intercepts by Group", x = "Group", y = "Effect (log-scale deviation)")

#    Emphasized on the importance of modeling group as a random effect to capture meaningful biological heterogeneity.

# Plot: Ocean random effects
ggplot(ocean_df_PFSAs, aes(x = reorder(ocean, effect), y = effect)) +
  geom_point(color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(title = "Random Intercepts by Ocean", x = "Ocean", y = "Effect (log-scale deviation)")

# Add a column to identify the random effect type
ocean_df_PFSAs$Type <- "Ocean"
group_df_PFSAs$Type <- "Genus"

# Rename the x-variable column to a common name, e.g. "Level"
ocean_df_PFSAs$Level <- ocean_df_PFSAs$ocean
group_df_PFSAs$Level <- group_df_PFSAs$group

# Combined Random Effects to plot together with two panels:

combined_df_PFSAs <- rbind(ocean_df_PFSAs[, c("Level", "effect", "Type")],
                           group_df_PFSAs[, c("Level", "effect", "Type")])


ggplot(combined_df_PFSAs, aes(x = reorder(Level, effect), y = effect)) +
  geom_point(color = "black", size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ Type, scales = "free_y") +
  labs(
    x = NULL,
    y = "Effect (log-scale deviation)") +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 10))

ggsave("random_effects_plot_PFSAs.pdf", width = 10, height = 6, dpi = 300)


  
