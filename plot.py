import numpy as np
import matplotlib.pyplot as plt




x, y = np.loadtxt("results\\particle_pos.txt").T

# x = x[len(x)/10:]
# y = y[len(y)/10:]
#
# x = x[::10]
# y = y[::10]

print("hi")

plt.plot(x,y)
plt.hold('on')

theta = np.linspace(0,2*np.pi,10000)

r = 1.0
x2 = r*np.cos(theta)
y2 = r*np.sin(theta)

plt.plot(x2, y2)

x, y = np.loadtxt("results\\planet_pos.txt").T

# x = x[len(x)/10:]
# y = y[len(y)/10:]
#
#
# x = x[::10]
# y = y[::10]
plt.plot(x,y)

ax = plt.gca()
ax.set_aspect('equal')

ax.set_xlim([-1.5,1.5])
ax.set_ylim([-1.5,1.5])


plt.show()
