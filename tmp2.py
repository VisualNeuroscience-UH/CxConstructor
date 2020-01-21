import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np

def func(input_currents,a,b,deltaI):
    # return a * np.exp(-b * x) + c
    return a*(input_currents - deltaI) + b*(input_currents - deltaI)**2
xdata = np.linspace(0, 4, 50)
y = func(xdata, 2.5, 1.3, 0.5)
np.random.seed(1729)
y_noise = 0.2 * np.random.normal(size=xdata.size)
ydata = y + y_noise
plt.plot(xdata, ydata, 'b-', label='data')
popt, pcov = curve_fit(func, xdata, ydata)
popt
plt.plot(xdata, func(xdata, *popt), 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.show()
