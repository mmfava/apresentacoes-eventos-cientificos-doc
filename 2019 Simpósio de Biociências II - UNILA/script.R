

## Trabajo accidentes Lonomia en Misiones

sexo=c("F"= 8, "M"=23)
chisq.test(sexo, simulate.p.value = T)
chisq.test(sexo, simulate.p.value = T)$residuals
(sexo/sum(sexo))*100

edad<-c(46, 22, 37, 3, 56, 14, 13, 34, 40, 5, 30, 12, 15, 11, 37, 15, 17, 9, 8, 9, 9, 6, 7, 17, 59, 2, 18, 59, 52, 12, 64)
boxplot(edad, ylab= "Edad")
summary(edad)
shapiro.test(edad)
sd(edad)
mean(edad)-sd(edad)

municipios=c("Bernardo de Irigoyen"=2, "Eldorado"=2, "Pozo Azul"=3, "Aristobulo del Valle"=2, 
             "San Pedro"=4, "San Vincente"=8, "San Antonio"=5, "25 de Mayo"=1, "Campo Ramón"=1, 
             "Colonia Aurora"=1, "El Soberbio"=2)
chisq.test(municipios, simulate.p.value = T)
chisq.test(municipios, simulate.p.value = T)$residuals


departamento=c("Gral Belgrano"=8, "Eldorado"= 2, "San Pedro"= 7, "Cainguas"= 2, "Guarani"= 9, "25 de Mayo"= 2, "Obera"=1)
chisq.test(departamento, simulate.p.value = T)
chisq.test(departamento, simulate.p.value = T)$residuals
sum(departamento)

ano=c("2014"= 7, "2015"= 12, "2016"= 6, "2017"=5, "2018"=1)
chisq.test(ano, simulate.p.value = T)
chisq.test(ano, simulate.p.value = T)$residuals
sum(ano)

mes=c("ene"= 6, "feb"=0, "mar"=0, "abr"=5, "may"=5, "jun"=0,
      "jul"=0, "ago"=0, "sep"=0, "oct"=0, "nov"=4, "dic"= 11)
chisq.test(mes, simulate.p.value = T)
chisq.test(mes, simulate.p.value = T)$residuals
sum(mes)

trime<-c(11, 5, 15)
chisq.test(trime)


tempo_atendimento=c("<1h"=3, "1-3h"= 1, "3-6h"= 1, "6-12h"=1, "12-24h"=3, ">24h"=22)
chisq.test(tempo_atendimento, simulate.p.value = T)
chisq.test(tempo_atendimento, simulate.p.value = T)$residuals
round((tempo_atendimento/(sum(tempo_atendimento))*100),2)

tempo_atend<-c(72, 0.96, 56, 24, 36, 24, 12, 48, 51, 3, 72, 12, 96, 24, 16, 24, 28, 72, 30, 10, 96, 159, 2, 24, 56, 28, 32, 28)
boxplot(tempo_atend)
summary(tempo_atend)
sd(tempo_atend)
shapiro.test(tempo_atend)

lugar=c("campo"=13 , "selva"=3 , "domicilio"=15 )
round(((lugar/sum(lugar))*100),2)
chisq.test(lugar, simulate.p.value = T)
chisq.test(lugar, simulate.p.value = T)$residuals

umbicación_picada=c("SD"=2, "multiples"= 1, "MMSS"= 20, "MMII"= 8)
round((umbicación_picada/sum(umbicación_picada))*100,2)
chisq.test(umbicación_picada, simulate.p.value = T)
chisq.test(umbicación_picada, simulate.p.value = T)$residuals

hemorragia=c("no"= 7, "sí"= 24)
chisq.test(hemorragia, simulate.p.value = T)
chisq.test(hemorragia, simulate.p.value = T)$residuals

classif<-c("leve"= 1, "moderado"= 28, "grave"= 2)
classif/sum(classif)*100
chisq.test(classif, simulate.p.value = T)
chisq.test(classif, simulate.p.value = T)$residuals

trat<-c("no"= 1, "5 ampollas"= 27, "8 ampollas"= 3)
chisq.test(trat, simulate.p.value = T)
chisq.test(trat, simulate.p.value = T)$residuals

ejemplar<-c("no"= 16, "foto"= 2, "traido"= 11, "capturado"= 2)
chisq.test(ejemplar, simulate.p.value = T)
chisq.test(ejemplar, simulate.p.value = T)$residuals


class_final<-c("descartado"= 8, "sospechoso"= 16, "confirmado"= 15)
round(((class_final/(sum(class_final)))*100),2)
chisq.test(class_final, simulate.p.value = T)


sintomas<-c(17, 6,	4,	2,	1,	1)
chisq.test(sintomas, simulate.p.value = T)
chisq.test(sintomas, simulate.p.value = T)$residuals



