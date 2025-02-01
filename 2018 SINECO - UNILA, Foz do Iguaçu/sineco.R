

## --- Script do SINECO --- ##

## ~~ Abrindo a base de dados ~~ ##
dados<-data.frame(dados)
names(dados)
attach(dados)

## ~~ Carregando os pacotes para as análises de dados ~~ ##
library(readr)
library(lme4)
library(ggplot2)
library(bbmle)
library(DHARMa)
library(AICcmodavg)
library(RVAideMemoire)
library(MASS)
library(INLA)
library(sjPlot)
library(visreg)

## ~~ Análises visuais (análises descritivas) ~~ ##


## Correlações 
plot(TA, TMC) 
cor.test(TA, TMC) ## 0.93

plot(delta_T1, delta_T2) 
cor.test(delta_T1, delta_T2)# -0.98

plot(PA, PMC)
cor.test(PA, PMC) # 0.40

plot(delta_P1, delta_P2)
cor.test(delta_P1, delta_P2) #-0.65

plot(delta_T1, delta_P1)
cor.test(delta_T1, delta_P1) #-0.60

plot(delta_T2, delta_P2)
cor.test(delta_T2, delta_P2) # -0.69

plot(delta_T1, TMC)
plot(delta_T2, TMC)

plot(S*100, AF*100)
plot(F,AF)
plot(S,F)

## ===== MODELOS ===== ##
names(dados)

## ~~~~~ POISSON ~~~~~ ##

# modelo nulo
PM0<-glm(n~1, data=dados, family=poisson(link=log)) 

# uso e ocupação da terra
PM1<-glm(n~S, data=dados, family=poisson(link=log)) 
PM2<-glm(n~AF, data=dados, family=poisson(link=log)) 
PM3<-glm(n~F, data=dados, family=poisson(link=log)) 
PM4<-glm(n~A, data=dados, family=poisson(link=log)) 

# clima
PM5<-glm(n~TMC, data=dados, family=poisson(link=log)) 
PM6<-glm(n~delta_T2, data=dados, family=poisson(link=log)) 
PM7<-glm(n~PMC, data=dados, family=poisson(link=log)) 
PM8<-glm(n~delta_P2, data=dados, family=poisson(link=log)) 

# Interação
PM9<-glm(n~S+AF, data=dados, family=poisson(link=log))
PM10<-glm(n~TMC+PMC, data=dados, family=poisson(link=log))
PM11<-glm(n~delta_T2+delta_P2, data=dados, family=poisson(link=log))
PM12<- glm(n~S+AF+TMC+PMC, data=dados, family=poisson(link=log))
PM13<- glm(n~S+AF+delta_T2+delta_P2, data=dados, family=poisson(link=log))
PM14<- glm(n~S+AF+TMC+PMC+delta_T2+delta_P2, data=dados, family=poisson(link=log))

## --> analise dos modelos poisson

poisson_models <- list(PM0, PM1, PM2, PM3, PM4, PM5, PM6, PM7, PM8, PM9, PM10, PM11, PM12, PM13, PM14)
names_PModels <- c("PM0","PM1","PM2","PM3", "PM4","PM5","PM6","PM7","PM8","PM9","PM10","PM11","PM12","PM13", "PM14")
print(aictab(cand.set = poisson_models, modnames = names_PModels, second.ord = F), digits = 3)
AIC(PM0, PM1, PM2, PM3, PM4, PM5, PM6, PM7, PM8, PM9, PM10, PM11, PM12, PM13, PM14)
BIC(PM0, PM1, PM2, PM3, PM4, PM5, PM6, PM7, PM8, PM9, PM10, PM11, PM12, PM13, PM14)
deviance(PM0)/df.residual(PM0) # 5.659134
deviance(PM1)/df.residual(PM1) # 5.655568
deviance(PM2)/df.residual(PM2) # 5.391663
deviance(PM3)/df.residual(PM3) # 5.662519
deviance(PM4)/df.residual(PM4) # 5.517128
deviance(PM5)/df.residual(PM5) # 5.517128
deviance(PM6)/df.residual(PM6) # 5.642927
deviance(PM7)/df.residual(PM7) # 5.549054
deviance(PM8)/df.residual(PM8) # 5.475722
deviance(PM9)/df.residual(PM0) # 5.377442
deviance(PM10)/df.residual(PM10) # 5.417975
deviance(PM11)/df.residual(PM11) # 5.427489
deviance(PM12)/df.residual(PM12) # 5.189895
deviance(PM13)/df.residual(PM13) # 5.189895
deviance(PM14)/df.residual(PM14) # 4.606114

## --> tabulando os modelos

# tabelando os modelos
tab_model(PM1, PM2, PM3, PM4, PM5, PM6, PM7, PM8, PM9, PM10, PM11, PM12, PM13, PM14, show.aic=T, show.ci=0.95, show.fstat=T)

# para visualizar os modelos 
visreg(PM14)
plotresid(PM14)
confint.merMod(PM14)
plot(PM14)
plot_model(PM14)

## ~~~~~ BINOMIAL NULO ~~~~~ ##

# modelo nulo
BN0<-glm.nb(n~1, data=dados) 

# uso e ocupação da terra
BN1<-glm.nb(n~S, data=dados) 
BN2<-glm.nb(n~AF, data=dados) 
BN3<-glm.nb(n~F, data=dados) 
BN4<-glm.nb(n~A, data=dados) 

# clima
BN5<-glm.nb(n~TMC, data=dados) 
BN6<-glm.nb(n~delta_T2, data=dados) 
BN7<-glm.nb(n~PMC, data=dados) 
BN8<-glm.nb(n~delta_P2, data=dados) 

# Interação
BN9<-glm.nb(n~S+AF, data=dados)
BN10<-glm.nb(n~TMC+PMC, data=dados)
BN11<-glm.nb(n~delta_T2+delta_P2, data=dados)
BN12<- glm.nb(n~S+AF+TMC+PMC, data=dados)
BN13<- glm.nb(n~S+AF+delta_T2+delta_P2, data=dados)
BN14<- glm.nb(n~S+AF+TMC+PMC+delta_T2+delta_P2, data=dados)

## --> analise dos modelos binomiais nulos

binomiais_nulos_models <- list(BN0, BN1, BN2, BN3, BN4, BN5, BN6, BN7, BN8, BN9, BN10, BN11, BN12, BN13)
names_BNModels <- c("BN0","BN1","BN2","BN3", "BN4","BN5","BN6","BN7","BN8","BN9","BN10","BN11","BN12","BN13")
AIC(BN0, BN1, BN2, BN3, BN4, BN5, BN6, BN7, BN8, BN9, BN10, BN11, BN12, BN13, BN14)
BIC(BN0, BN1, BN2, BN3, BN4, BN5, BN6, BN7, BN8, BN9, BN10, BN11, BN12, BN13, BN14)

deviance(BN0)/df.residual(BN0) # 0.4752219
deviance(BN1)/df.residual(BN1) # 0.4757579
deviance(BN2)/df.residual(BN2) # 0.4799047
deviance(BN3)/df.residual(BN3) # 0.4756434
deviance(BN4)/df.residual(BN4) # 0.4777303
deviance(BN5)/df.residual(BN5) # 0.4777303
deviance(BN6)/df.residual(BN6) # 0.4758713
deviance(BN7)/df.residual(BN7) # 0.4777714
deviance(BN8)/df.residual(BN8) # 0.4788611
deviance(BN9)/df.residual(BN0) # 0.4796222
deviance(BN10)/df.residual(BN10) # 0.4798938
deviance(BN11)/df.residual(BN11) # 0.4798938
deviance(BN12)/df.residual(BN12) # 0.484541
deviance(BN13)/df.residual(BN13) # 0.484541
deviance(BN14)/df.residual(BN14) # 0.4921662


## --> tabulando os modelos

# tabelando os modelos
tab_model(BN1, BN2, BN3, BN4, BN5, BN6, BN7, BN8, BN9, BN10, BN11, BN12, BN13, BN14, show.aic=T, show.ci=0.95, show.fstat=T)

# para visualizar os modelos 
visreg(BN14)
plotresid(BN14)
plot(BN14)
plot_model(BN14)

## ~~~~~ BINOMIAL ~~~~~ ##

# modelo nulo
BM0<-glm(pres~1, data=dados, family=binomial(link=logit)) 

# uso e ocupação da terra
BM1<-glm(pres~S, data=dados, family=binomial(link=logit)) 
BM2<-glm(pres~AF, data=dados, family=binomial(link=logit)) 
BM3<-glm(pres~F, data=dados, family=binomial(link=logit)) 
BM4<-glm(pres~A, data=dados, family=binomial(link=logit)) 

# clima
BM5<-glm(pres~TMC, data=dados, family=binomial(link=logit)) 
BM6<-glm(pres~delta_T2, data=dados, family=binomial(link=logit)) 
BM7<-glm(pres~PMC, data=dados, family=binomial(link=logit)) 
BM8<-glm(pres~delta_P2, data=dados, family=binomial(link=logit)) 

# Interação
BM9<-glm(pres~S+AF, data=dados, family=binomial(link=logit))
BM10<-glm(pres~TMC+PMC, data=dados, family=binomial(link=logit))
BM11<-glm(pres~delta_T2+delta_P2, data=dados, family=binomial(link=logit))
BM12<- glm(pres~S+AF+TMC+PMC, data=dados, family=binomial(link=logit))
BM13<- glm(pres~S+AF+delta_T2+delta_P2, data=dados, family=binomial(link=logit))
BM14<- glm(pres~S+AF+TMC+PMC+delta_T2+delta_P2, data=dados, family=binomial(link=logit))


## --> analise dos modelos binomiais

binomial_models <- list(BM0, BM1, BM2, BM3, BM4, BM6, BM8, BM9, BM11, BM13)
names_BModels <- c("BM0","BM1","BM2","BM3", "BM4","BM6","BM8","BM9","BM11","BM13")

print(aictab(cand.set = binomial_models, modnames = names_BModels, second.ord = F), digits = 3)
AIC(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14)
BIC(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14)


deviance(BM0)/df.residual(BM0) # 1.037375
deviance(BM1)/df.residual(BM1) # 1.033415
deviance(BM2)/df.residual(BM2) # 1.024869
deviance(BM3)/df.residual(BM3) # 1.038084
deviance(BM4)/df.residual(BM4) # 1.019636
deviance(BM5)/df.residual(BM5) # 1.011154
deviance(BM6)/df.residual(BM6) # 1.036352
deviance(BM7)/df.residual(BM7) # 1.018803
deviance(BM8)/df.residual(BM8) # 1.0321
deviance(BM9)/df.residual(BM9) # 1.012128
deviance(BM10)/df.residual(BM10) # 0.9960885
deviance(BM11)/df.residual(BM11) # 1.032742
deviance(BM12)/df.residual(BM12) # 0.9887564
deviance(BM13)/df.residual(BM13) # 1.00759
deviance(BM14)/df.residual(BM14) # 0.9426938


## --> tabulando os modelos

# tabelando os modelos
tab_model(BM1, BM2, BM3, BM4, BM6,  BM8, BM9, BM11, BM13, show.aic=T, show.ci=0.95, show.fstat=T)

# para visualizar os modelos 
visreg(BM13)
plotresid(BM13)
plot(BM13)
plot_model(BM13)


##### ------------- MODELOS GENERALIZADOS E MIXTOS ------------- #####

## ~~~~~ Binomial nulo ~~~~~ ##

## Modelo nulo 
BMM0<-glmer(pres~1+(1|Cidade), data=dados, family=binomial(link="logit")) 

# uso e ocupação da terra
BMM1<-glmer(pres~S+(1|Cidade), data=dados, family=binomial(link="logit")) 
BMM2<-glmer(pres~AF+(1|Cidade), data=dados, family=binomial(link="logit")) 
BMM3<-glmer(pres~F+(1|Cidade), data=dados, family=binomial(link="logit")) 
BMM4<-glmer(pres~A+(1|Cidade), data=dados, family=binomial(link="logit")) 

# clima
BMM5<-glmer(pres~TMC+(1|Cidade), data=dados, family=binomial(link="logit")) 
BMM6<-glmer(pres~delta_T2+(1|Cidade), data=dados, family=binomial(link="logit")) 
BMM7<-glmer(pres~PMC+(1|Cidade), data=dados, family=binomial(link="logit")) 
BMM8<-glmer(pres~delta_P2+(1|Cidade), data=dados, family=binomial(link="logit")) 

# Interação
BMM9<-glmer(pres~S+AF+(1|Cidade), data=dados, family=binomial(link="logit"))
BMM10<-glmer(pres~TMC+PMC+(1|Cidade), data=dados, family=binomial(link="logit"))
BMM11<-glmer(pres~delta_T2+delta_P2+(1|Cidade), data=dados, family=binomial(link="logit"))
BMM12<- glmer(pres~S+AF+TMC+PMC+(1|Cidade), data=dados, family=binomial(link="logit"))
BMM13<- glmer(pres~S+AF+delta_T2+delta_P2+(1|Cidade), data=dados, family=binomial(link="logit"))
BMM14<- glmer(pres~S+AF+TMC+PMC+delta_T2+delta_P2+(1|Cidade), data=dados, family=binomial(link="logit"))

## --> analise dos modelos binomiais mixtos
binomial_mixmodels <- list(BMM0, BMM1, BMM2,BMM3,BMM4,BMM5,BMM6,BMM7,BMM8,BMM9,BMM10,BMM11,BMM12,BMM13,BMM14)
names_BMModels <- c("BMM0", "BMM1", "BMM2","BMM3","BMM4","BMM5","BMM6","BMM7","BMM8","BMM9","BMM10","BMM11","BMM12","BMM13","BMM14")
print(aictab(cand.set = binomial_models, modnames = names_BModels, second.ord = F), digits = 3)
AIC(BMM0, BMM1, BMM2,BMM3,BMM4,BMM5,BMM6,BMM7,BMM8,BMM9,BMM10,BMM11,BMM12,BMM13,BMM14)
BIC(BMM0, BMM1, BMM2,BMM3,BMM4,BMM5,BMM6,BMM7,BMM8,BMM9,BMM10,BMM11,BMM12,BMM13,BMM14)

## --> analise dos modelos binomiais
deviance(BMM0)/df.residual(BMM0) # 1.038248
deviance(BMM1)/df.residual(BMM1) # 1.034284
deviance(BMM2)/df.residual(BMM2) # 0.9935517
deviance(BMM3)/df.residual(BMM3) # 1.038958
deviance(BMM4)/df.residual(BMM4) # 1.020494
deviance(BMM5)/df.residual(BMM5) # 0.9343354
deviance(BMM6)/df.residual(BMM6) # 1.024194
deviance(BMM7)/df.residual(BMM7) # 0.9410279
deviance(BMM8)/df.residual(BMM8) # 1.0321
deviance(BMM9)/df.residual(BMM9) # 1.031989
deviance(BMM10)/df.residual(BMM10) # 0.9216873
deviance(BMM11)/df.residual(BMM11) # 1.033612
deviance(BMM12)/df.residual(BMM12) # 0.935018
deviance(BMM13)/df.residual(BMM13) # 0.970423

## ~~~~~ Poisson ~~~~~ ##

## Modelo nulo 
PMM0<-glmer(n~1+(1|Cidade), data=dados, family=poisson(link="log")) 

# uso e ocupação da terra
PMM1<-glmer(pres~S+(1|Cidade), data=dados, family=poisson(link="log")) 
PMM2<-glmer(pres~AF+(1|Cidade), data=dados, family=poisson(link="log")) 
PMM3<-glmer(pres~F+(1|Cidade), data=dados, family=poisson(link="log")) 
PMM4<-glmer(pres~A+(1|Cidade), data=dados, family=poisson(link="log")) 

# clima
PMM5<-glmer(pres~TMC+(1|Cidade), data=dados, family=poisson(link="log")) 
PMM6<-glmer(pres~delta_T2+(1|Cidade), data=dados, family=poisson(link="log")) 
PMM7<-glmer(pres~PMC+(1|Cidade), data=dados, family=poisson(link="log")) 
PMM8<-glmer(pres~delta_P2+(1|Cidade), data=dados, family=poisson(link="log")) 

# Interação
PMM9<-glmer(pres~S+AF+(1|Cidade), data=dados, family=poisson(link="log"))
PMM10<-glmer(pres~TMC+PMC+(1|Cidade), data=dados, family=poisson(link="log"))
PMM11<-glmer(pres~delta_T2+delta_P2+(1|Cidade), data=dados, family=poisson(link="log"))
PMM12<- glmer(pres~S+AF+TMC+PMC+(1|Cidade), data=dados, family=poisson(link="log"))
PMM13<- glmer(pres~S+AF+delta_T2+delta_P2+(1|Cidade), data=dados, family=poisson(link="log"))
PMM14<- glmer(pres~S+AF+TMC+PMC+delta_T2+delta_P2+(1|Cidade), data=dados, family=poisson(link="log"))

## --> analise dos modelos binomiais mixtos
poisson_mixmodels <- list(PMM0, PMM1, PMM2,PMM3,PMM4,PMM5,PMM6,PMM7,PMM8,PMM9,PMM10,PMM11,PMM12,PMM13,PMM14)
names_PMModels <- c("PMM0", "PMM1", "PMM2","PMM3","PMM4","PMM5","PMM6","PMM7","PMM8","PMM9","PMM10","PMM11","PMM12","PMM13","PMM14")
print(aictab(cand.set = binomial_models, modnames = names_BModels, second.ord = F), digits = 3)
AIC(PMM0, PMM1, PMM2,PMM3,PMM4,PMM5,PMM6,PMM7,PMM8,PMM9,PMM10,PMM11,PMM12,PMM13,PMM14)
BIC(PMM0, PMM1, PMM2,PMM3,PMM4,PMM5,PMM6,PMM7,PMM8,PMM9,PMM10,PMM11,PMM12,PMM13,PMM14)

## --> analise dos modelos poisson
deviance(PMM0)/df.residual(PMM0) # 0.04641932
deviance(PMM1)/df.residual(PMM1) # 0.6570633
deviance(PMM2)/df.residual(PMM2) # 0.6501871
deviance(PMM3)/df.residual(PMM3) # 0.6606211
deviance(PMM4)/df.residual(PMM4) # 0.6458002
deviance(PMM5)/df.residual(PMM5) # 0.6399873
deviance(PMM6)/df.residual(PMM6) # 0.6592586
deviance(PMM7)/df.residual(PMM7) # 0.6461909
deviance(PMM8)/df.residual(PMM8) # 0.6558784
deviance(PMM9)/df.residual(PMM9) # 0.6399679
deviance(PMM10)/df.residual(PMM10) # 0.6279971
deviance(PMM11)/df.residual(PMM11) # 0.6562426
deviance(PMM12)/df.residual(PMM12) # 0.6215113
deviance(PMM13)/df.residual(PMM13) # 0.6361582

## ~~~~~ Poisson ~~~~~ ##

## Modelo nulo 
BNM0<-glmer.nb(n~1+(1|Cidade), data=dados) 

# uso e ocupação da terra
BNM1<-glmer.nb(pres~S+(1|Cidade), data=dados) 
BNM2<-glmer.nb(pres~AF+(1|Cidade), data=dados) 
BNM3<-glmer.nb(pres~F+(1|Cidade), data=dados) 
BNM4<-glmer.nb(pres~A+(1|Cidade), data=dados) 

# clima
BNM5<-glmer.nb(pres~TMC+(1|Cidade), data=dados) 
BNM6<-glmer.nb(pres~delta_T2+(1|Cidade), data=dados) 
BNM7<-glmer.nb(pres~PMC+(1|Cidade), data=dados) 
BNM8<-glmer.nb(pres~delta_P2+(1|Cidade), data=dados) 

# Interação
BNM9<-glmer.nb(pres~S+AF+(1|Cidade), data=dados)
BNM10<-glmer.nb(pres~TMC+PMC+(1|Cidade), data=dados)
BNM11<-glmer.nb(pres~delta_T2+delta_P2+(1|Cidade), data=dados)
BNM12<- glmer.nb(pres~S+AF+TMC+PMC+(1|Cidade), data=dados)
BNM13<- glmer.nb(pres~S+AF+delta_T2+delta_P2+(1|Cidade), data=dados)
BNM14<- glmer.nb(pres~S+AF+TMC+PMC+delta_T2+delta_P2+(1|Cidade), data=dados)

## --> analise dos modelos binomiais mixtos
poisson_mixmodels <- list(BNM0, BNM1, BNM2,BNM3,BNM4,BNM5,BNM6,BNM7,BNM8,BNM9,BNM10,BNM11,BNM12,BNM13,BNM14)
names_BNModels <- c("BNM0", "BNM1", "BNM2","BNM3","BNM4","BNM5","BNM6","BNM7","BNM8","BNM9","BNM10","BNM11","BNM12","BNM13","BNM14")
print(aictab(cand.set = binomial_models, modnames = names_BModels, second.ord = F), digits = 3)
AIC(BNM0, BNM1, BNM2,BNM3,BNM4,BNM5,BNM6,BNM7,BNM8,BNM9,BNM10,BNM11,BNM12,BNM13,BNM14)
BIC(BNM0, BNM1, BNM2,BNM3,BNM4,BNM5,BNM6,BNM7,BNM8,BNM9,BNM10,BNM11,BNM12,BNM13,BNM14)

## --> analise dos modelos poisson
deviance(BNM0)/df.residual(BNM0) # 0.04625026
deviance(BNM1)/df.residual(BNM1) # 0.6576054
deviance(BNM2)/df.residual(BNM2) # 0.6507233
deviance(BNM3)/df.residual(BNM3) # 0.6611661
deviance(BNM4)/df.residual(BNM4) # 0.6463319
deviance(BNM5)/df.residual(BNM5) # 0.640515
deviance(BNM6)/df.residual(BNM6) # 0.6598025
deviance(BNM7)/df.residual(BNM7) # 0.6467237
deviance(BNM8)/df.residual(BNM8) # 0.6564194
deviance(BNM9)/df.residual(BNM9) # 0.6404956
deviance(BNM10)/df.residual(BNM10) # 0.6285148
deviance(BNM11)/df.residual(BNM11) # 0.6567843
deviance(BNM12)/df.residual(BNM12) # 0.6220239
deviance(BNM13)/df.residual(BNM13) # 0.6366835

#### --------------------------------------------------- ####

binomial_models <- list(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14)
names_BModels <- c("BM0","BM1","BM2","BM3", "BM4","BM5", "BM6", "BM7", "BM8","BM9","BM10", "BM11", "BM12", "BM13", "BM14")

print(aictab(cand.set = binomial_models, modnames = names_BModels, second.ord = F), digits = 3)
AIC(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14)
BIC(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14)
deviance(BM14)/df.residual(BM14) # 0.9426938
summary(BM14)
tab_model(BM14)
logLik(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14)
logLik(BM0)
logLik(BM14)
logLik(BM13)


plot_model(BM14, axis.labels=c("Agricola|Floresta", "??P", "??T", "PMC", "Silvicultura", "TMC"))
plot_model(BM14, type="est")
plot_model(BM14, type="std2")
plot_model(BM14, type="resid")
plot_model(BM14, type="diag")

visreg(BM14)

aic_bm<-list(AIC(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14))
aic_bm

bic_bm<-list(BIC(BM0, BM1, BM2, BM3, BM4, BM5, BM6, BM7, BM8, BM9, BM10, BM11, BM12, BM13,BM14))
bic_bm

dev_bm<-list(deviance(BM0)/df.residual(BM0),
             deviance(BM1)/df.residual(BM1),
             deviance(BM2)/df.residual(BM2),
             deviance(BM3)/df.residual(BM3),
             deviance(BM4)/df.residual(BM4),
             deviance(BM5)/df.residual(BM5),
             deviance(BM6)/df.residual(BM6),
             deviance(BM7)/df.residual(BM7),
             deviance(BM8)/df.residual(BM8),
             deviance(BM9)/df.residual(BM9),
             deviance(BM10)/df.residual(BM10),
             deviance(BM11)/df.residual(BM11),
             deviance(BM12)/df.residual(BM12),
             deviance(BM13)/df.residual(BM13),
             deviance(BM14)/df.residual(BM14))

dev_bm

bb<-data.frame(aic_bm, bic_bm)
bb

fixef(BM14)
profile(BM14)

plot(predict(BM14), residuals(BM14))
summary(BM14)

min(S)
min(AF)
min(TMC)
min(PMC)
min(delta_T2)
min(delta_P2)

mean(S)
mean(AF)
mean(TMC)
mean(PMC)
mean(delta_T2)
mean(delta_P2)

max(S)
max(AF)
max(TMC)
max(PMC)
max(delta_T2)
max(delta_P2)

n
