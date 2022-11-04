import numpy as np
import random as rd
import pandas as pd
from scipy import optimize
from matplotlib import pyplot as plt


def generate_data(gantry):
    x, data, count = [], pd.DataFrame(columns=['gantry', 'dx']), 20
    def randsin(x):  return np.sin(2*np.pi / 360 * x) + rd.gauss(0., 1.)
    for i in range(count):
        x.extend(gantry)

    data.gantry = x
    data.dx = data.gantry.apply(randsin)
    
    return data


def fit_sin(t, y):
    tt, yy = np.array(t), np.array(y)
    
    def sinfunc(t, A, w, p, c):  return A * np.sin(w*t + p) + c

    guess_freq = 1 / 360
    guess_amp = np.std(yy) * 2.**0.5
    guess_offset = np.mean(yy)
    guess = np.array([guess_amp, 2.*np.pi*guess_freq, 0., guess_offset])

    popt, pcov = optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*np.pi)
    fitfunc = lambda t: A * np.sin(w*t + p) + c

    print('Amplitude:', A, '\nPeriod:', 1/f, '\nPhase:', p, '\nOffset:', c)
    return fitfunc


rd.seed(999)
gantry = np.linspace(0., 360., 8, endpoint=False)
x_axis = np.linspace(0, 315, 1000)
data = generate_data(gantry)
means = [data.loc[data.gantry == gtr].dx.mean() for gtr in gantry]
fit = fit_sin(gantry, means)

data.plot('gantry', 'dx', kind='scatter')
plt.scatter(gantry, means, c='r')
plt.plot(x_axis, fit(x_axis), c='r')
plt.show()