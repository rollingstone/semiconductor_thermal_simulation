__author__ = 'kamal'



from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import random

count = []

for i in range(1, 4):
    for j in range(3, 6):
        for k in range(15,19):
            count.append((i, j, k, random.random()))

data = np.array(count)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(data[:,0], data[:,1], data[:,3])
plt.show()

