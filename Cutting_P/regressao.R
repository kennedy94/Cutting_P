
data = read.csv(file.choose())
#plot(data)

modelo = lm(DESEM ~ TP + NG + MUT + RST + AS + CRS + TER,data)
modelo = lm(SNR ~ TP + NG + MUT + RST + AS + CRS + TER,data)
summary(anova(modelo))
#aov(DESEM ~ TP + NG + MUT + RST + AS + CRS + TER,data)

plot(data$Instance,data$LP, type = "+")
plot(data$LB, type = "o", lty=1, pch=1, xlab = "Instance",ylab = "Value")
lines(data$IP, type = "o", lty=2, pch=3)
lines(data$LP, type = "o", lty=5, pch = 2)
legend("topright", legend = c(" ", " "))