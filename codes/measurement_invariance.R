# required packages
install.packages(c("haven", "lavaan", "semPlot", "knitr", "dplyr", "tidyr"))

library(haven) # to read spss files
library(lavaan) # to run cfa()
library(semPlot) # to run semPaths()
library(knitr)
library(dplyr)
library(tidyr)

# reading data
PISA2015 <- read_sav("datawithoutoutlier.sav")

# data manipulation
PISA2015[, c("OZYET_1", "OZYET_2", "OZYET_3", "OZYET_4", "OZYET_5", "OZYET_6", "OZYET_7", "OZYET_8",
             "MOT_1", "MOT_2", "MOT_3", "MOT_4")] <- lapply(PISA2015[, c("OZYET_1", "OZYET_2", "OZYET_3",
                                                                         "OZYET_4", "OZYET_5", "OZYET_6",
                                                                         "OZYET_7", "OZYET_8", "MOT_1",
                                                                         "MOT_2", "MOT_3", "MOT_4")], ordered)

PISA2015$BOLGE <- as.factor(PISA2015$BOLGE)
PISA2015$CINSIYET <- as.factor(PISA2015$CINSIYET)

# model fitting
model <- 'OZYET =~ NA*OZYET_1 + OZYET_2 + OZYET_3 + OZYET_4 + 1*OZYET_5 + OZYET_6 + OZYET_7 + OZYET_8
          MOT =~ NA*MOT_1 + 2*MOT_2 + MOT_3 + MOT_4'

fit <- cfa(model, data = PISA2015, estimator = "ULS")
           
summary(fit, fit.measures = TRUE, standardized = TRUE)

# goodness of fit indexes
fitMeasures(fit, c("rmsea", "cfi", "tli", "srmr", "gfi"))

# subsetting for cfa         
PISA2015_KIZ <- subset(PISA2015, CINSIYET == 1)
PISA2015_ERKEK <- subset(PISA2015, CINSIYET == 2)

fit_kiz <- cfa(model, PISA2015_KIZ, estimator = "ULS")
fit_kiz
fit_erkek <-cfa(model, PISA2015_ERKEK, estimator = "ULS")
fit_erkek

# goodness of fit indexes
fitMeasures(fit_kiz , c("rmsea", "cfi", "tli", "srmr", "gfi"))
fitMeasures(fit_erkek , c("rmsea", "cfi", "tli", "srmr", "gfi"))

# confirmatory factor analysis visualization
semPaths(fit, "std", rotation = 2, layout = "tree2", nCharNodes = 0, 
         sizeLat = 8, sizeLat2 = 6, mar = c(2, 6, 2, 4), edge.color = "black")

# factor loadings
options(knitr.kable.NA = '')
parameterEstimates(fit, standardized = TRUE) |> 
  filter(op == "=~") |> 
  mutate(stars = ifelse(pvalue < .001, "***", 
                        ifelse(pvalue < .01, "**", 
                               ifelse(pvalue < .05, "*", "")))) |>
  select('Latent Factor' = lhs, 
         Indicator = rhs, 
         B = est, 
         SE = se, Z = z, 
         Beta = std.all, 
         sig = stars) |> 
  kable(digits = 3, format = "pandoc", caption = "Table 1: Factor Loadings")


# grouping by CINSIYET
## configural model
fit_configural_cinsiyet <- cfa(model, data = PISA2015, estimator = "ULS", group = "CINSIYET")
fit_configural_cinsiyet

## weak model
fit_weak_cinsiyet <- cfa(model, data = PISA2015, estimator ="ULS", group = "CINSIYET", group.equal = "loadings")
fit_weak_cinsiyet

## configural model versus weak model
anova(fit_weak_cinsiyet, fit_configural_cinsiyet)

## strong model
fit_strong_cinsiyet <- cfa(model, data = PISA2015, estimator = "ULS", group = "CINSIYET",
                           group.equal = c("loadings", "intercepts"))
fit_strong_cinsiyet

## weak model vs strong model
anova(fit_strong_cinsiyet, fit_weak_cinsiyet)

## strict model
fit_strict_cinsiyet <- cfa(model, data = PISA2015, estimator = "ULS", group = "CINSIYET",
                           group.equal = c("loadings", "intercepts", "residuals"))
fit_strict_cinsiyet

## goodness of fit indexes
fit.stats_cinsiyet <- rbind(fitMeasures(fit_configural_cinsiyet,
                                        c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")),
                            fitMeasures(fit_weak_cinsiyet,
                                        c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")),
                            fitMeasures(fit_strong_cinsiyet,
                                        c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")),
                            fitMeasures(fit_strict_cinsiyet,
                                        c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")))
fit.stats_cinsiyet

# chi-squared difference test
lavTestLRT(fit_configural_cinsiyet, fit_weak_cinsiyet, fit_strong_cinsiyet, fit_strict_cinsiyet)




# grouping by BOLGE
## confirmatory factor analysis
## BOLGE == 1
PISA2015_B1 <- subset(PISA2015, BOLGE == 1)
fit_bolge1 <- cfa(model, PISA2015_B1, estimator = "ULS")
fit_bolge1
fitMeasures(fit_bolge1, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 2
PISA2015_B2 <- subset(PISA2015, BOLGE == 2)
fit_bolge2 <- cfa(model, PISA2015_B2, estimator = "ULS")
fit_bolge2
fitMeasures(fit_bolge2, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 3
PISA2015_B3 <- subset(PISA2015, BOLGE == 3)
fit_bolge3 <- cfa(model, PISA2015_B3, estimator = "ULS")
fit_bolge3
fitMeasures(fit_bolge3, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 4
PISA2015_B4 <- subset(PISA2015, BOLGE == 4)
fit_bolge4 <- cfa(model, PISA2015_B4, estimator = "ULS")
fit_bolge4
fitMeasures(fit_bolge4, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 5
PISA2015_B5 <- subset(PISA2015, BOLGE == 5)
fit_bolge5 <- cfa(model, PISA2015_B5, estimator = "ULS")
fit_bolge5
fitMeasures(fit_bolge5, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 6
PISA2015_B6 <- subset(PISA2015, BOLGE == 6)
fit_bolge6 <- cfa(model, PISA2015_B6, estimator = "ULS")
fit_bolge6
fitMeasures(fit_bolge6, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 7
PISA2015_B7 <- subset(PISA2015, BOLGE == 7)
fit_bolge7 <- cfa(model, PISA2015_B7, estimator = "ULS")
fit_bolge7
fitMeasures(fit_bolge7, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 8
PISA2015_B8 <- subset(PISA2015, BOLGE == 8)
fit_bolge8 <- cfa(model, PISA2015_B8, estimator = "ULS")
fit_bolge8
fitMeasures(fit_bolge8, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 9
PISA2015_B9 <- subset(PISA2015, BOLGE == 9)
fit_bolge9 <- cfa(model, PISA2015_B9, estimator = "ULS")
fit_bolge9
fitMeasures(fit_bolge9, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 10
PISA2015_B10 <- subset(PISA2015, BOLGE == 10)
fit_bolge10 <- cfa(model, PISA2015_B10, estimator = "ULS")
fit_bolge10
fitMeasures(fit_bolge10, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 11
PISA2015_B11 <- subset(PISA2015, BOLGE == 11)
fit_bolge11 <- cfa(model, PISA2015_B11, estimator = "ULS")
fit_bolge11
fitMeasures(fit_bolge11, c("rmsea", "cfi", "tli", "srmr", "gfi"))

## BOLGE == 12
PISA2015_B12 <- subset(PISA2015, BOLGE == 12)
fit_bolge12 <- cfa(model, PISA2015_B12, estimator = "ULS")
fit_bolge12
fitMeasures(fit_bolge12, c("rmsea", "cfi", "tli", "srmr", "gfi"))

# measurement invariance analysis will be continued with the new data set
PISA2015_Bolge <-  subset(PISA2015, BOLGE == 1 | BOLGE == 3 | BOLGE == 5 | BOLGE == 6 | BOLGE == 11 | BOLGE == 12)               

## configural model
fit_configural_bolge <- cfa(model, data = PISA2015_Bolge, estimator = "ULS", group = "BOLGE")
fit_configural_bolge

## weak model
fit_weak_bolge <- cfa(model, data = PISA2015_Bolge, estimator = "ULS", group = "BOLGE", group.equal = "loadings")
fit_weak_bolge

## configural model versus weak model
anova(fit_weak_bolge, fit_configural_bolge)

## strong model
fit_strong_bolge <- cfa(model, data = PISA2015_Bolge, estimator = "ULS", group = "BOLGE", group.equal = c("loadings", "intercepts"))
fit_strong_bolge

## weak model vs strong model
anova(fit_strong_bolge, fit_weak_bolge)

## strict model
fit_strict_bolge <- cfa(model, data = PISA2015_Bolge, estimator = "ULS", group = "BOLGE",
                        group.equal = c("loadings", "intercepts", "residuals"))
fit_strict_bolge

## goodness of fit indexes
fit.stats_bolge <- rbind(fitMeasures(fit_configural_bolge, 
                                     c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")),
                         fitMeasures(fit_weak_bolge,
                                     c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")),
                         fitMeasures(fit_strong_bolge,
                                     c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")),
                         fitMeasures(fit_strict_bolge,
                                     c("chisq", "rmsea", "cfi", "tli", "rni", "rfi", "ifi", "srmr", "gfi")))
fit.stats_bolge

## chi-squared difference test
lavTestLRT(fit_configural_bolge, fit_weak_bolge, fit_strong_bolge, fit_strict_bolge)
