import matplotlib.pyplot as plt
import numpy as np


def f1(x, a, b):
    return a + b * x


def f2(x, a, b, c):
    return a + b * x + c * x * x


coef = []
plt.grid(True)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.xticks([-3, -2, -1, 0, 1, 2])
plt.yticks([-2, -1, 0, 1, 2, 3, 4, 5, 6, 7])
with open('C:/Users/pavel/Desktop/DA/NM_3_3/test.txt') as f:
    for line in f:
        s = line.rstrip()
        s = s.split(" ")
        for val in s:
            convert = float(val)
            coef.append(convert)


x = [-3, -2, -1, 0, 1, 2]
y = [0.04979, 0.13534, 0.36788, 1, 2.7183, 7.3891]

xlist = np.linspace(-3, 2, num=1000)
ylist = f1(xlist, coef[0], coef[1])
flist = f2(xlist, coef[2], coef[3], coef[4])

plt.scatter(x, y, marker='o', color="red")

plt.xlabel('x - axis')
plt.ylabel('y - axis')

plt.plot(xlist, ylist)
plt.plot(xlist, flist)
plt.show()
