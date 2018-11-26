
data = read.csv(file.choose())
library(phia)
library(tree)
plot(data)

modelo = lm(SNR2 ~ TP + NG + MUT + TER, data)
modelo = lm(SNR ~ TP + NG + MUT + TER,data)
modelo = lm(SNRT ~ NG + RST + AS,data)

summary(aov(modelo))
summary(anova(modelo))
#aov(DESEM ~ TP + NG + MUT + RST + AS + CRS + TER,data)

#plot
layout(matrix(c(1,2,3,4,5,6,7,7), 2, 4, byrow = TRUE))
boxplot(DESEM ~ TP, data, xlab="TP")
boxplot(DESEM ~ NG, data, xlab="NG")
boxplot(DESEM ~ RST, data, xlab="RST")
boxplot(DESEM ~ AS, data, xlab="AS")
boxplot(DESEM ~ CRS, data, xlab="CRS")
boxplot(DESEM ~ TER, data, xlab="TER")
boxplot(DESEM ~ MUT, data, xlab="MUT")

