from mpl_chord_diagram import chord_diagram
import numpy as np
import matplotlib.pyplot as plt

flux_data = np.array([
    [0, 0.001, 0.030, 0.001, 0.001, 0.001, 0.001, 0.001],
    [0.001, 0, 0.036, 0.001, 0.034, 0.038, 0.001, 0.028],
    [0.030, 0.036, 0, 0.001, 0.001, 0.036, 0.001, 0.001],
    [0.001, 0.001, 0.001, 0, 0.001, 0.036, 0.001, 0.001],
    [0.001, 0.034, 0.001, 0.001, 0, 0.033, 0.001, 0.001],
    [0.001, 0.038, 0.036, 0.036, 0.033, 0, 0.001, 0.001],
    [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0, 0.001],
    [0.001, 0.028, 0.001, 0.001, 0.001, 0.001, 0.001, 0],
])
names = ["PZ","P3","P4","POZ","O1","O2","PO7","PO8"]
filtered_flux_data = flux_data * (flux_data > 0.001)

fig,ax = plt.subplots(figsize=(4,3.5),dpi=100,facecolor="w")
chord_diagram = chord_diagram(mat=filtered_flux_data,names=names,alpha=.8,ax=ax, use_gradient=True)
for text in ax.texts:
    text.set_fontsize(5) 
plt.tight_layout()
plt.show()