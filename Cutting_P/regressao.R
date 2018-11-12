
data = read.csv(file.choose())
#plot(data)

modelo = lm(DESEM ~ TP + NG + MUT + RST + AS + CRS + TER,data)
modelo = lm(SNR ~ TP + NG + MUT + RST + AS + CRS + TER,data)
anova(modelo)
summary(anova(modelo))
#aov(DESEM ~ TP + NG + MUT + RST + AS + CRS + TER,data)
interaction.plot(x.factor = data$MUT,
                 response = data$DESEM,
                 trace.factor = data$CRS,
                 fun = mean,
                 col=c("black","red","green"),  ### Colors for levels of trace var.
                 pch=c(19, 17, 15),             ### Symbols for levels of trace var.
                 fixed=TRUE,                    ### Order by factor order in data
                 leg.bty = "o")