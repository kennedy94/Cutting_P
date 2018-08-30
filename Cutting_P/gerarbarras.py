import pandas as pd
import io
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties

import sys

df = pd.read_csv(sys.argv[1], header=None, names=["Task", "Start", "Finish", "Resource"] )
df["Diff"] = df.Finish - df.Start


color = {"Type 1":"0.25", "Type 2":"0.5", "Type 3":"0", "Type 0":"1"} #, "JOB3":"blue", "JOB4":"red"}


blackpatch	= mpatches.Patch(color='0.25', label='Bar')
graypatch	= mpatches.Patch(color='0.5', label='Leftover')
whitepatch	= mpatches.Patch(color='0', label='Waste')
#bluepatch	= mpatches.Patch(color='blue', label='JOB4')
#red_patch	= mpatches.Patch(color='red', label='JOB5')

fig,ax=plt.subplots(figsize=(6,3))


labels=[]


for i, task in enumerate(df.groupby("Task")):
	labels.append(task[0])
	for r in task[1].groupby("Resource"):
		data = r[1][["Start", "Diff"]]
		ax.broken_barh(data.values, (i-0.4,0.8), color=color[r[0]] )



ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels) 
#ax.set_xlabel("Time")

fontP = FontProperties()
fontP.set_size('small')
#plt.legend(handles=[blackpatch,graypatch,greenpatch,bluepatch, red_patch],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.legend(handles=[blackpatch,graypatch,whitepatch],bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.tight_layout()       
#plt.show()

nome  = sys.argv[1]

if nome.endswith('.solu'):
	nome = nome[:-5]
	nome = nome + '_gantt.png'
if nome.endswith('.barras'):
	nome = nome[:-7]
	nome = nome + '_barras.png'
	


plt.savefig(nome, bbox_inches='tight')