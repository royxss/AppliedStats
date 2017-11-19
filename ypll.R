rm(list=ls())
seedVal = 17869
options(warn=-1)
options(scipen=999)

library(ggplot2)
library(ggcorrplot)
library(dummies)
library(rms)
library(lmtest)

measures <- read.csv2("additional_measures_cleaned.csv", header = TRUE, sep = ',')
ypll <- read.csv2("ypll.csv", header = TRUE, sep = ',')

head(measures)
head(ypll)
df <- merge(x = ypll, y = measures, by.x = c("FIPS"), by.y = c("FIPS"))
str(df)
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
df <- df[(df$County.x)!='',]

apply(df, 2, function(x) length(which(x == "" | is.na(x) | x == "NA")))
# HIV rate is an important predictor. but 22% values are missing
# so it is difficult to impute
# we can bucket it into categories
HIVRate <- df[!is.na(df$HIV.rate),'HIV.rate']
range(HIVRate)
quantile(HIVRate)

# we need to convert HiveRate to %HIV multiplying it by 100 more as binning requires int not num
percentHIV <- HIVRate/df$Population*100*100
range(percentHIV)
quantile(percentHIV)
length(percentHIV)
length(df$FIPS)
#df$HIV <- cut(df$HIV.rate, breaks = c(-Inf, 57, 98, 194, Inf), 
#              labels = c("VeryLow","Low","High","VeryHigh"))

df$HIV <- cut(percentHIV, breaks = c(-Inf, 12, 34, 94, Inf), 
              labels = c("VeryLow","Low","High","VeryHigh"))

df[is.na(df$HIV),'HIV'] <- 'NotAvail'
summary(df)
names(df)
str(df)

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


# Outliers

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
modelForm <- createModelFormula(yVar = yVar, xVars = xVars, includeIntercept = TRUE)
modelForm
model <- lm(modelForm, data = df)
summary(model)
plot(model)


#Test for constant variance
summary(lm(abs(residuals(model)) ~ fitted(model)))

#It turns out that the absolute residuals are not predicted very well by Yˆi
#Hence, we conclude that there does not seem to be a problem with the constant variance assumption.

#Test for Normality
hist(residuals(model))
boxplot(residuals(model))

#From the plots it seems like the residuals are normally distributed with mean of 0

#Shapiro.Wilks Normality Test
shapiro.test(residuals(model))
#The null hypothesis is that the residuals have a normal distribution. The p-value of the test
#statistic is large in this example. It thus follows that the null hypothesis is not rejected.

# to check for serial correlation
#Durbin-Watson Test
dwtest(model,data=df)
# Th null hypothesis is there is no correlation among errors. The p value indicates that there is no evidence of correlated errors, 

# Test stepwise
model_step <- step(model, direction="both")
summary(model_step)

plot(model_step)

#Test for constant variance
summary(lm(abs(residuals(model_step)) ~ fitted(model_step)))

#It turns out that the absolute residuals are not predicted very well by Yˆi
#Hence, we conclude that there does not seem to be a problem with the constant variance assumption.

#Test for Normality
hist(residuals(model_step))
boxplot(residuals(model_step))

#From the plots it seems like the residuals are normally distributed with mean of 0

#Shapiro.Wilks Normality Test
shapiro.test(residuals(model_step))
#The null hypothesis is that the residuals have a normal distribution. The p-value of the test
#statistic is large in this example. It thus follows that the null hypothesis is not rejected.


# Fit model with transformation
logTransformModelForm <- createModelFormula(yVar = 'log(YPLL.Rate)', xVars = xVars, includeIntercept = TRUE)
logTransformModelForm

logTransform <- lm(logTransformModelForm, data = df)
#summary(glm(modelForm,data=df))
summary(logTransform)
plot(logTransform)

#Test for constant variance
summary(lm(abs(residuals(logTransform)) ~ fitted(logTransform)))

#It turns out that the absolute residuals are not predicted very well by Yˆi
#Hence, we conclude that there does not seem to be a problem with the constant variance assumption.

#Test for Normality
hist(residuals(logTransform))
boxplot(residuals(logTransform))

#From the plots it seems like the residuals are normally distributed with mean of 0

#Shapiro.Wilks Normality Test
shapiro.test(residuals(logTransform))
#The null hypothesis is that the residuals have a normal distribution. The p-value of the test
#statistic is large in this example. It thus follows that the null hypothesis is not rejected.

# to check for serial correlation 
#Durbin-Watson Test
#install.packages('lmtest')
library(lmtest)
dwtest(logTransform,data=df)
# The null hypothesis is there is no correlation among errors. The p value indicates that there is no evidence of correlated errors, 


#Leverage

#hatvalues(model)

hv <- as.data.frame(hatvalues(model))
hvmean <-mean(hatvalues(model))
hv$warn <- ifelse(hv[, 'hatvalues(model)']>3*hvmean, 'x3',
                  ifelse(hv[, 'hatvalues(model)']>2*hvmean, 'x2', '-' ))

plot(hatvalues(model), type = "h")

x3<- which(hv$warn == "x3")
x2 <- which(hv$warn %in%c("x2", "x3"))

leveragex3subsetdf <- df[-x3,]
hv3model <- lm(modelForm, data = leveragex3subsetdf)
summary(hv3model)
plot(hv3model)


leveragex2subsetdf <- df[-x2,]
hv2model <- lm(modelForm, data = leveragex2subsetdf)
summary(hv2model)
plot(hv2model)



#Inluence TestCooks D
#cooks.distance(model)

cooksD <- as.data.frame(cooks.distance(model))
cooksmean <-mean(cooks.distance(model))
cooksD$warn <- ifelse(cooksD[, 'cooks.distance(model)']>4*cooksmean, 'x4','-' )

plot(cooksD$`cooks.distance(model)`,pch="*", cex=2, main="Influential Obs by Cooks distance")
abline(h = 4*cooksmean, col="red")
x4<- which(cooksD$warn == "x4")

inluencex4subsetdf <- df[-x4,]
cooksmodel <- lm(modelForm, data = inluencex4subsetdf)
summary(cooksmodel)
plot(cooksmodel)

