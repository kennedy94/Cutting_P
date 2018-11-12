library(AlgDesign)
levels.design = c(4,3,3,3)
f.design <- gen.factorial(levels = levels.design,
                          varNames = c("Popu","Gen", "Mut", "Restart", "NAlea", "Crossover","tax_elitismo"),
                          factors = "all")

fract.design <- optFederov(
  data=f.design,
  nTrials=10,
  #approximate=TRUE,
  criterion = "D")

print(fract.design)
print(f.design)