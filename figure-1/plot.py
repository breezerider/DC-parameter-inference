import matplotlib.pyplot as plt
import numpy as np


gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))

cmap = plt.get_cmap('RdYlGn')
cmap_r = cmap.reversed()

fig = plt.figure(frameon=False)
fig.set_size_inches(8,1)
ax = plt.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)

ax.imshow(gradient, aspect='auto', cmap=cmap_r)

fig.savefig(f'colorbar.png', dpi=96.0)
plt.close(fig)

names = ['false-negatives_costfcn','false-negatives_oracle','false-positives_costfcn','false-positives_oracle']

for n, f in [('homo', 'Homo-Reaction-Model'), ('sbd', 'Single-Birth-Death-Model')]:
  a = np.loadtxt(f'comparison_{f}.csv',delimiter=',')
  a[:,0:2] /= 9900.0
  a[:,2:4] /= 10000.0
  for i in range(4):
    b = 1.0 - a[:,i]
    b = b.reshape((100,100))
    
    fig = plt.figure(frameon=False)
    fig.set_size_inches(8,8)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    ax.imshow(b.T, interpolation='bicubic', cmap=cmap_r, aspect='auto')
    
    fig.savefig(f'{names[i]}_{n}.png', dpi=192.0)
    plt.close(fig)
    
    #plt.show()
  
