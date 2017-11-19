setwd("C:\\Users\\SROY\\Desktop\\Courses\\Applied Stats\\Project")
rm(list=ls())
seedVal = 17869
options(warn=-1)
options(scipen=999)

# Libraries
library(ggplot2)
library(ggcorrplot)
library(dummies)

measures <- read.csv2("additional_measures_cleaned.csv", header = TRUE, sep = ',')
ypll <- read.csv2("ypll.csv", header = TRUE, sep = ',')

str(measures)
str(ypll)

df <- merge(x = ypll, y = measures, by.x = c("FIPS"), by.y = c("FIPS"))

apply(df, 2, function(x) length(which(x == "" | is.na(x) | x == "NA")))

# remove records that are unreliable
df <- df[df$Unreliable != 'x',names(df) != 'Unreliable']
# Remove cumulative records
df <- df[df$County.x != '',]
# Remove missing YPLL rate as it does not add value and 
# doesn't make sense to impute it as it is the target variable
df <- df[!is.na(df$YPLL.Rate),]
# Remove missing Rural, X.free.lunch and X..chile.Illiteracy as they are few
# removing 1 rural also removes free lunch and child illiteracy 
df <- df[!df$Rural == '',]
df <- df[!df$X..child.Illiteracy == '',]
df <- df[!is.na(df$X..Free.lunch),]

# HIV rate is an important predictor. but 22% values are missing
# so it is difficult to impute
# we can bucket it into categories
HIVRate <- df[!is.na(df$HIV.rate),'HIV.rate']
range(HIVRate)
quantile(HIVRate)

df$HIV <- cut(df$HIV.rate, breaks = c(-Inf, 57, 98, 194, Inf), 
              labels = c("VeryLow","Low","High","VeryHigh"))
df[is.na(df$HIV),'HIV'] <- 'NotAvail'

# Remove irrelevant columns
names(df)
excludeList <- c("FIPS", "State.x", "State.y", "County.x", "County.y", "HIV.rate")
includeList <- names(df[!names(df) %in% c(excludeList)])
df <- df[,includeList]

# Change factors to numeric
df$X..18 <- as.numeric(as.character(df$X..18))
df$X65.and.over <- as.numeric(as.character(df$X65.and.over))
df$African.American <- as.numeric(as.character(df$African.American))
df$Female <- as.numeric(as.character(df$Female))
df$Rural <- as.numeric(as.character(df$Rural))
df$X..child.Illiteracy <- as.numeric(as.character(df$X..child.Illiteracy))
str(df)
summary(df)

# select numeric columns
numVars <- names(which(sapply(df, is.numeric)))
# correlation matrix
corr <- round(cor(df[,numVars]), 1)
ggcorrplot(corr, hc.order = TRUE, 
           type = "lower", 
           lab = TRUE, 
           lab_size = 3, 
           method="circle", 
           colors = c("tomato2", "white", "springgreen3"), 
           title="Correlation matrix", 
           ggtheme=theme_bw)

# Histogram
ggplot(df, aes(YPLL.Rate)) +  
  geom_histogram(fill="tomato3") +  
  labs(title="Histogram of YPLL Rate")
# not highly skewed. No need to scale. We can verify it later.

# Scatterplot
ggplot(df, aes(x=Population, y=YPLL.Rate)) + 
  geom_point() + 
  geom_smooth(method="loess", se=F) + scale_x_log10() +
  labs(x="Population scaled to log 10")

# Dummy HIV variable
df <- cbind(df, dummy(df$HIV))
df <- df[,!names(df) %in% "HIV"]

# Assign x y vars
yVar <- "YPLL.Rate"
xVars <- names(df[!names(df) %in% c(yVar)])

createModelFormula <- function(yVar, xVars, includeIntercept = TRUE){
  if(includeIntercept){
    modelForm <- as.formula(paste(yVar, "~", paste(xVars, collapse = '+ ')))
  } else {
    modelForm <- as.formula(paste(yVar, "~", paste(xVars, collapse = '+ '), -1))
  }
  return(modelForm)
}

# Fit linear model with all variables
modelForm <- createModelFormula(yVar = yVar, xVars = xVars, includeIntercept = FALSE)
modelForm
lm.1 <- lm(modelForm, data = df)
#summary(lm.1)
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.1)


# Remove outliers using studentized residual
stud_resid <- rstudent(lm.1)
hist(stud_resid, main = "Studentized Residual Distribution")
# Exclude studentized residual > 3
rmv_stud_resid <- stud_resid[abs(stud_resid) > 2]
length(rmv_stud_resid)


# Exclude outliers
df_exc_out <- df[!row.names(df) %in% names(rmv_stud_resid),]
lm.2 <- lm(modelForm, data = df_exc_out)
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.2)

# Reset rownames
rownames(df) <- seq(length=nrow(df))
rownames(df_exc_out) <- seq(length=nrow(df_exc_out))

# Remove influential observations using Cook's D
cook_simple = cooks.distance(lm.1)
cook_exc_outliers = cooks.distance(lm.2)
# identify Cook's D values > 4/(n-k-1) 
cutoff <- 4/((nrow(df_exc_out)-length(lm.2$coefficients)-2))
par(mfrow = c(1, 2))
plot(as.integer(rownames(df)), cook_simple, main = "Including Outliers", 
     ylab = "Cook's Distance", ylim = c(0, max(cook_simple)))
abline(h = 0.005)
plot(as.integer(rownames(df_exc_out)), cook_exc_outliers, main = "Excluding Outliers", 
     ylab = "", ylim = c(0, max(cook_simple)))
abline(h = 0.005)

# Find Cook on top of data without outliers
rmv_cook_inf <- cook_exc_outliers[cook_exc_outliers > cutoff]
length(rmv_cook_inf)

# Exclude Influence
df_exc_inf <- df_exc_out[!row.names(df_exc_out) %in% names(rmv_cook_inf),]
lm.3 <- lm(modelForm, data = df_exc_inf)
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.3)


# Scale values to reduce Heteroskedasticity
str(df_exc_inf)
df_exc_scalex <- df_exc_inf
# We can scale only x variables to ln (Population and median.household.income) and verify
df_exc_scalex$Population <- log(df_exc_scalex$Population, exp(1))
df_exc_scalex$median.household.income <- log(df_exc_scalex$median.household.income, exp(1))
lm.4 <- lm(modelForm, data = df_exc_scalex)
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.4)

# We now can scale only y to ln (YPLL.Rate) and verify
df_exc_scaley <- df_exc_inf
df_exc_scaley$YPLL.Rate <- log(df_exc_scaley$YPLL.Rate, exp(1))
lm.5 <- lm(modelForm, data = df_exc_scaley)
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.5) 


# We now can scale both x and y to ln and verify
df_exc_scalexy <- df_exc_inf
df_exc_scalexy$Population <- log(df_exc_scalexy$Population, exp(1))
df_exc_scalexy$median.household.income <- log(df_exc_scalexy$median.household.income, exp(1))
df_exc_scalexy$YPLL.Rate <- log(df_exc_scalexy$YPLL.Rate, exp(1))
lm.6 <- lm(modelForm, data = df_exc_scalexy)
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.6) 

# Test stepwise
lm.7 <- step(lm.5, direction="both")
# Residual analysis:
par(mfrow = c(2, 2))
plot(lm.7)

# We now know model 5 is the best. But does step wise improve model 5?
# Compare lm 7 with lm 5 as it is built on top of it.
# Reduction in the residual sum of squares are statistically significant or not)
anova(lm.5,lm.7, test="Chisq")

# We choose model 7 as the best model as it gives the same performance with less features
summary(lm.5)
summary(lm.7)


