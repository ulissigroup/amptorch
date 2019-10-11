import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 1, 100).reshape(-1, 1)
fc = 0.5*(1+np.cos(np.pi*x))
rs = 0
t = 80
g2 = np.exp(-t*(x**2))*fc
fig, ax = plt.subplots(figsize=(14.15, 10))
plt.plot(x, g2 ,color='k', linewidth=3)
plt.xlim(left=0, right=1)
plt.ylim(bottom=0, top=1)
plt.xlabel(r'$r/R_{c}$', fontsize=25)
plt.ylabel(r'$f_{c}$', fontsize=25)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.show()
# plt.savefig('cos_cutoff.png',dpi=300)
