import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pylab import *


data = loadtxt("Baseball.txt")
mtofeet = 3.2808
t = data[:,0]
x = data[:,1]#*mtofeet
y = data[:,2]#*mtofeet
z = data[:,3]#*mtofeet
vx = data[:,4]
vy = data[:,5]
vz = data[:,6]
#print(data)
#print(t)
#print(x)
#print(v)
fig = plt.figure(1)

#plt.plot(t,x, color="red",linestyle="-.")
#plt.plot(t,y, color="black",linestyle="-.")
#plt.plot(t,z, color="black",linestyle="-.")


plt.plot(x,y, color="black",linestyle="-.")
plt.plot(x,z, color="blue",linestyle="-.")

plt.legend(["y","z"])
#plt.legend(["M2","E2","M4","E4"])
#plt.legend(["M2","E2","M3","E3","M4","E4"])
plt.title("Baseball: Slider, right-handed")
plt.xlabel("x/m")
plt.ylabel("distance/m")
fig.savefig('Slider.png', bbox_inches='tight')
#plt.ylabel("")
plt.show()

