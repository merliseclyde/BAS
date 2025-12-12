library(BAS)
data(UScrime, package="MASS")
#UScrime[,-2] = log(UScrime[,-2])

set.seed(42)
grow.time = system.time(crime.grow <-  bas.lm(log(y) ~ log(M) + So + log(Ed) + log(Po1) + log(Po2)
                      + log(LF) + log(M.F) + log(Pop) + log(NW) +
                        log(U1) + log(U2) + log(GDP) + log(Ineq) + log(Prob)+
                        log(Time), 
                      data=UScrime, n.models=2^15, prior="BIC", 
                      method = "MCMC", burnin.iteration = 5000, 
                      MCMC.iterations = 500000, GROW = TRUE))
print(c(crime.grow$n.models, grow.time))
set.seed(42)
old.time = system.time(crime.old <-  bas.lm(log(y) ~ log(M) + So + log(Ed) + log(Po1) + log(Po2)
                                              + log(LF) + log(M.F) + log(Pop) + log(NW) +
                                                log(U1) + log(U2) + log(GDP) + log(Ineq) + log(Prob)+
                                                log(Time), 
                                              data=UScrime, n.models=2^15, prior="BIC", 
                                              method = "MCMC", burnin.iteration = 5000, 
                                              MCMC.iterations = 500000, GROW = FALSE, n.models.init = 2^15))

print(c(crime.old$n.models, old.time))


data(tecator, package="FuncNN")
# Extract data and target
X <- tecator$absorp.fdata$data[1:172,]
fat <- tecator$y$Fat[1:172] 
data = data.frame(fat, X)

subsamp = seq(1,100, by=4)
p = length(subsamp)
b.it = 20000
mc.it = 200000
n.models = 2^20

set.seed(42)
old.time = system.time(tec.old <-  bas.lm(fat ~ ., data=data[, c(subsamp,101)],
                                            n.models=n.models, prior="BIC", 
                                            method = "MCMC", burnin.iteration = b.it, 
                                            MCMC.iterations = mc.it, GROW = FALSE, n.models.init = n.models))

print(c(tec.old$n.models, old.time))

set.seed(42)
grow.time = system.time(tec.grow <-  bas.lm(fat ~ ., data=data[, c(subsamp,101)],
                                          n.models=n.models, prior="BIC", 
                                          method = "MCMC", burnin.iteration = b.it, 
                                          MCMC.iterations = mc.it, GROW = TRUE))

print(c(tec.grow$n.models, grow.time))
