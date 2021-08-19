import matplotlib.pyplot as plt
import numpy as np
import sys

def readf(f):
    x = []
    y = []
    for line in f:
        a = float(line.split()[0])      
        x.append(a)
        y.append(float(line.split()[1]))
    return x, y


f = open ('data2.txt', 'r')
lda = np.arange (600, 800, 10)
thick = np.arange(40e-9, 60e-9, 1e-9)
z = np.zeros((len(lda), len(thick)))
i = 0
for line in f:
  for j in range(len(thick)):
    z[i][j] = float(line.split()[j])
  i = i + 1

print (z.shape)
L, H = np.meshgrid(lda, thick)
fig, ax = plt.subplots(1, 1, figsize = (20, 20))
s = ax.contourf(lda, thick, z, cmap='hot')
ax.set_ylabel ("thickness, m")
ax.set_xlabel ("wavelength, nm")
fig.colorbar(s)
fig.savefig("map2.png")
plt.close()