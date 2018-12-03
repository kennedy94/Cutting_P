# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 10:49:40 2018

@author: Kennedy
"""

import matplotlib.pyplot as plt
#import numpy as np
import pandas as pd

dados = pd.read_csv("medias.csv")
dados.head()
plot = dados.plot(kind = "line", grid = True,
           color = ['gray','black'], style = '-')
plt.xscale('symlog', linthreshy=0.015)
#plt.yticks(np.arange(0, 4000, step=500))
#plt.show()
plt.savefig('tempo_ga_cplex_log.png', dpi=250)





dados = pd.read_csv("tempos_ga_cplex.csv")
plot = dados.plot(kind = "line", grid = True,
           color = ['gray','black'], style = '.-')
plt.yscale('symlog', linthreshy=0.015)
#plt.yticks(np.arange(0, 4000, step=500))
plot.set_xlabel('Instance')
plot.set_ylabel('time (s)')
plt.savefig('tempo_ga_cplex_log.png', dpi=250)
plt.show()



dados = pd.read_csv("tempos_ga_cplex.csv")
plot = dados.plot(kind = "line", grid = True,
           color = ['gray','black'], style = '.-')
#plt.yscale('symlog', linthreshy=0.015)
#plt.yticks(np.arange(0, 4000, step=500))
plot.set_xlabel('Instance')
plot.set_ylabel('time (s)')
plt.savefig('tempo_ga_cplex.png', dpi=250)
plt.show()



dados = pd.read_csv("comparacao_GA_CPLEX.csv")
dados.head()
plot = dados.plot(kind = "line", grid = True,
           color = ['gray','black'], style = '.-')

plot.set_xlabel('Instance')
plot.set_ylabel('LBD')

plt.savefig('GAvsCPLEX.png', dpi=250)


dados = pd.read_csv("LB_IP_LP.csv")

plot = dados.plot(kind = "line", grid = True,
           color = [ 'red', 'gray','black'],
           style = ['.-','+-','*-'])
plot.set_xlabel('Instance')
plot.set_ylabel('Objective function value')

plt.savefig('lowerboundCP.png', dpi=250)


#import seaborn
#seaborn.pairplot(dados, vars=['LB'],
#                 kind='reg', hue='IP')  

dados = pd.read_csv("data.csv")

dados.head()

#dados.boxplot(column = 'DESEM', 
#              by = ['TP','NG','MUT','RST','AS','CRS','TER'], grid = False)


fig, axes = plt.subplots(nrows=2, ncols=4, figsize=(15, 10))

dados.boxplot(column = 'DESEM', 
              by = 'TP', grid = False,
              ax=axes[0,0])


dados.boxplot(column = 'DESEM', 
              by = 'NG', grid = False,
              ax=axes[0,1])

dados.boxplot(column = 'DESEM', 
              by = 'MUT', grid = False,
              ax=axes[0,2])

dados.boxplot(column = 'DESEM', 
              by = 'RST', grid = False,
              ax=axes[0,3])

dados.boxplot(column = 'DESEM', 
              by = 'AS', grid = False,
              ax=axes[1,0])

dados.boxplot(column = 'DESEM', 
              by = 'CRS', grid = False,
              ax=axes[1,1])

dados.boxplot(column = 'DESEM', 
              by = 'TER', grid = False,
              ax=axes[1,2])

axes[1, 3].axis('off')

plt.savefig('teste.png', dpi=250)



