library(devtools)
# devtools::install_github("betsybersson/BAS")
# library(BAS)
library(dplyr)

devtools::load_all()

# sim data
n=50
Y = rgamma(n,2,2)+runif(n)/2
X1 = rnorm(n,2,5)
plot(X1,Y)
hist(Y)
df = cbind(Y,X1) %>% as.data.frame()


# BAS model
out =  bas.glm(Y~ ., data=df, family = Gamma(link = "log"))
out$mle

## compare to base R model
out_freq = glm(Y~.,data=df,family=Gamma(link="log"))
coef(out_freq,dispersion = 0.1)

## better example
library(faraway)
data(wafer)
wafer_glm <- glm(formula = resist ~ .,
                 family  = Gamma(link = "log"),
                 data    = wafer)

wafer_bas = bas.glm(resist~ ., data=wafer, family = Gamma(link = "log"))
