library(AlgDesign)
levels.design = c(3,2,3,2,2,2,2)
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