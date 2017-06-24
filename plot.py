import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Malth.txt")

print(data)

t = data[:,0]
L = data[:,1]
Y = data[:,2]
y = data[:,3]
Lind = data[:,4]
Yind = data[:,5]
yind = data[:,6]




plt.plot(t,Lind)
plt.plot(t,Yind)
plt.plot(t,yind)
plt.legend(["L = population","Y = GDP","y = GDP/capita"])
plt.xlabel("time t")
plt.title("Economic and population model based on Malthusian principles, Index t=0")
plt.ylabel("Index")
plt.savefig('figs/Malthus.png', bbox_inches='tight')
plt.show()
