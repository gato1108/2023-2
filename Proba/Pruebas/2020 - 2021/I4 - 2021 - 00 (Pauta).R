############################
## I4 - 2021 - 00 (Pauta) ##
############################

library(TeachingDemos)
library(rio)
library(dplyr)


##################
## Pregunta 1.1 ##
##################

## H0: p = p0 vs Ha: p > p0
p0 = 1/3
n = 60
X = 25
hat.p = X/n

c(z.test(x = hat.p, mu = p0, stdev = sqrt(p0*(1-p0)/n), alternative = "greater")$statistic,
z.test(x = hat.p, mu = p0, stdev = sqrt(p0*(1-p0)/n), alternative = "greater")$p.value)

c(prop.test(x = X, n = n, p = p0, correct = F, alternative = "greater")$statistic,
prop.test(x = X, n = n, p = p0, correct = F, alternative = "greater")$p.value)

Z0 = (hat.p-p0)/sqrt(p0*(1-p0)/n)
c(Z0,1-pnorm(Z0))

## +0.5 por Z0
## +0.3 por valor-p
## +0.2 por conclusión: NO

##################
## Pregunta 1.2 ##
##################

## H0: mu = mu0 vs Ha: mu < mu0
mu0 = 5.0
n = 60-25
hat.mu = 4.7
s = 1.2
T0 = (hat.mu-mu0)/(s/sqrt(n))
c(T0, pt(T0,df = n-1))

## +0.5 por T0
## +0.3 por valor-p
## +0.2 por conclusión: NO

################
## Pregunta 2 ##
################

## Ejemplo
p = 0.39
omega = 0.05
confianza = 0.91
alpha = 1-confianza
n = trunc((qnorm(1-alpha/2)*sqrt(p*(1-p))/omega)^2)+1
n

## +1.0 por n +/-2

##################
## Pregunta 3.1 ##
##################

## H0: p_c = p_h vs Ha: p_c != p_h
c(prop.test(x = c(32, 21), n = c(72, 62), correct = F, alternative = "two.sided")$statistic,
prop.test(x = c(32, 21), n = c(72, 62), correct = F, alternative = "two.sided")$p.value)

p = (32+21)/(72+62)
Z0 = c(32/72-21/62)/(sqrt(p*(1-p))*sqrt(1/72+1/62))
c(Z0,2*(1-pnorm(abs(Z0))))

c(z.test(x = 32/72-21/62, stdev = (sqrt(p*(1-p))*sqrt(1/72+1/62)))$statistic,
z.test(x = 32/72-21/62, stdev = (sqrt(p*(1-p))*sqrt(1/72+1/62)))$p.value)

## +0.5 por estadistico
## +0.3 por valor-p
## +0.2 por conclusión: NO

##################
## Pregunta 3.2 ##
##################

## H0: mu_c = mu_h vs Ha: mu_c != mu_h
F0 = 1.2^2/0.9^2
## F0>1
2*(1-pf(F0, df1 = 32-1, df2 = 21-1))
sp = sqrt(((32-1)*1.2^2+(21-1)*0.9^2)/(32+21-2))
T0 = c(4.9-5.4)/(sp*sqrt(1/32+1/21))
c(T0, 2*(1-pt(abs(T0), df = 32+21-2)))

## +0.5 por estadistico
## +0.3 por valor-p
## +0.2 por conclusión: NO

T0 = c(4.9-5.4)/(sqrt(1.2^2/32+0.9^2/21))
nu = (1.2^2/32+0.9^2/21)^2/((1.2^2/32)^2/(32-1)+(0.9^2/21)^2/(21-1))
c(T0, 2*(1-pt(abs(T0), df = nu)))

## +0.3 por estadistico incorrecto
## +0.3 por valor-p
## +0.2 por conclusión: NO

################
## Pregunta 4 ##
################

## H0: mu = 8 vs Ha: mu < 8

Base = as.data.frame(rio::import("Salario_I4.xlsx"))
X = Base$Experiencia
c(t.test(x = X, mu = 8, alternative = "less")$statistic,
t.test(x = X, mu = 8, alternative = "less")$p.value)

## +0.5 por estadistico
## +0.3 por valor-p
## +0.2 por conclusión: NO

################
## Pregunta 5 ##
################

Base = as.data.frame(rio::import("Salario_I4.xlsx"))
X = dplyr::filter(Base, Genero == "MASCULINO")$Ingreso
Y = dplyr::filter(Base, Genero == "FEMENINO")$Ingreso
var.test(x = X, y = Y, alternative = "two.sided")$p.value
c(t.test(x = X, y = Y, alternative = "greater", var.equal = T)$statistic,
t.test(x = X, y = Y, alternative = "greater", var.equal = T)$p.value)

## +0.5 por estadistico
## +0.3 por valor-p
## +0.2 por conclusión: NO

c(t.test(x = X, y = Y, alternative = "greater", var.equal = F)$statistic,
t.test(x = X, y = Y, alternative = "greater", var.equal = F)$p.value)

## +0.3 por estadistico incorrecto
## +0.3 por valor-p
## +0.2 por conclusión: NO

################
## Pregunta 6 ##
################

Base = as.data.frame(rio::import("Salario_I4.xlsx"))
modelo1 = lm(Ingreso ~ Experiencia, data = Base)
modelo2 = lm(Ingreso ~ Experiencia + I(Experiencia^2) + Genero, data = Base)

## Función anova() permite comparar modelos anidados (gerarquicos) 
## y evaluar el aporte de las nuevas variables ingresadas como un conjunto
anova(modelo1,modelo2)

## +0.4 por estadistico
## +0.4 por valor-p
## +0.2 por conclusión: NO

################
## Pregunta 7 ##
################

## hat.nu ~ Normal(nu, nu/sqrt(n))
hat.nu = 1/3.4
n = 133
alpha = 1-0.90
LI = hat.nu - qnorm(1-alpha/2)*1/sqrt(n)
LS = hat.nu + qnorm(1-alpha/2)*hat.nu/sqrt(n)
cbind(LI,LS)

## 1.0 por IC sobre nu (Modelo Exponencial)

hat.mu = 3.4
n = 133
alpha = 1-0.90
LI = hat.mu - qnorm(1-alpha/2)*sqrt(hat.mu)/sqrt(n)
Ls = hat.mu + qnorm(1-alpha/2)*sqrt(hat.mu)/sqrt(n)
cbind(LI,LS)

## 0.8 por IC sobre mu (Modelo Poisson) [Incorrecto]

################
## Pregunta 8 ##
################

## H0: 1/3 cada profesor vs Ha: distinto a 1/3 al menos uno
p = c(1/3,1/3,1/3)
O = c(18,10,8)
chisq.test(x = O, p = p)
E = c(12,12,12)
chisq.test(x = O, p = E, rescale.p = T)

## +0.4 por estadistico
## +0.4 por valor-p
## +0.2 por conclusión: SI

##################
## Pregunta 9.1 ##
##################

library(readxl)
library(fitdistrplus)
Base = as.data.frame(read_excel("Palta_2015_2020.xlsx"))
X = Base$precio
par = fitdist(data = X, distr = "lnorm", method = "mle")$estimate
ks.test(X, "plnorm", meanlog = par[1], sdlog = par[2])$p.value

## +0.8 por valor-p
## +0.2 por conclusión: NO

##################
## Pregunta 9.2 ##
##################

library(readxl)
library(fitdistrplus)
Base = as.data.frame(read_excel("Palta_2015_2020.xlsx"))
X = Base$precio
par = fitdist(data = X, distr = "weibull", method = "mle")$estimate
ks.test(X, "pweibull", shape = par[1], scale = par[2])$p.value

## +0.8 por valor-p
## +0.2 por conclusión: SI

##################
## Pregunta 9.3 ##
##################

library(readxl)
library(fitdistrplus)
Base = as.data.frame(read_excel("Palta_2015_2020.xlsx"))
X = Base$precio
par = fitdist(data = X, distr = "gamma", method = "mle")$estimate
ks.test(X, "pgamma", shape = par[1], rate = par[2])$p.value

## +0.8 por valor-p
## +0.2 por conclusión: NO

#################
## Pregunta 10 ##
#################

library(readxl)
library(dplyr)
Base = as.data.frame(read_excel("Palta_2015_2020.xlsx"))
X = filter(Base, tipo == "organic")$precio
xp = sort(X)
N = length(X)
p = 1:N/(N+1)
par = lm(log(xp) ~ log(-log(1-p)))$coef
eta = exp(par[1])
beta = 1/par[2]
cbind(eta, beta)

## +1.0 por valores correctos

par = lm((xp) ~ log(-log(1-p)))$coef
eta = exp(par[1])
beta = 1/par[2]
cbind(eta, beta)

## +0.5 por valores correctos al olvidar aplicar log(xp)
