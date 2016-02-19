__author__ = 'Justin Finkle'
__email__ = 'jfinkle@u.northwestern.edu'

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib as mpl
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['font.sans-serif'] = 'Arial'
# mpl.rcParams['xtick.major.pad'] = 12
# mpl.rcParams['ytick.major.pad'] = 12


def biochem_sim_poly(x):
    base_x = np.arange(36)
    base_y = np.array([1, 1, 1, 1, 1, 1, 2, 4, 6, 8, 10, 12, 14, 15, 16, 17, 18, 18, 17.5, 17, 16, 15, 13, 11, 8, 7,
                       6.5, 6.25, 6, 6, 6, 6, 6, 6, 6, 6])
    poly = np.poly1d(np.polyfit(base_x, base_y, 10))
    y = poly(x)
    start_min = np.min(y[:len(y)/2])
    start_min_index = np.where(y==start_min)[0][0]
    stop_min = np.min((y[len(y)/2:]))
    stop_min_index = np.where(y==stop_min)[0][0]
    y[:start_min_index] = start_min
    y[stop_min_index:] = stop_min
    return y

def shift_values(x, y, noise, shift):

    points_to_prepend = np.sum(x < shift)
    shifted_noise = np.append(noise[:points_to_prepend], noise)
    shifted_y = np.append(np.array([y[0]]*points_to_prepend), y)+shifted_noise

    return shifted_y[:len(y)]

if __name__ == '__main__':
    cutoffs = range(0, 11, 5)[::-1]
    np.random.seed(10)
    n_points = 100
    noise_factor = 1
    gauss1 = np.random.normal(0, 1, n_points)*noise_factor
    gauss2 = np.random.normal(0, 1, n_points)*noise_factor
    a = np.linspace(0, 36, n_points)
    b = biochem_sim_poly(a)
    plot_b = stats.zscore(b+gauss1, ddof=1)
    line_width = 2
    f, axarr = plt.subplots(2, len(cutoffs), figsize=(15, 10))
    tick_size = 18

    for col, cutoff in enumerate(cutoffs):

        b2 = stats.zscore(shift_values(a, b, gauss2, cutoff), ddof=1)

        # Fit linear model
        slope, intercept, r_value, p_value, std_err = stats.linregress(b, b2)

        axarr[0, col].plot(a, plot_b, lw=line_width, label='Gene 1', c='b')
        axarr[0, col].plot(a, b2, lw=line_width, label='Gene 2', c='r')
        axarr[0, col].set_xlim([np.min(a), np.max(a)])
        axarr[0, col].set_title('Lag Order: %i' % (10-cutoff), fontsize=tick_size, weight='bold')
        axarr[0, col].tick_params(axis='both', labelsize=tick_size)
        if col != 0:
            axarr[0, col].get_xaxis().set_visible(False)
            axarr[0, col].get_yaxis().set_visible(False)
        else:
            axarr[0, col].legend(loc='best')
            axarr[0, col].locator_params(nbins=4)

        axarr[1, col].plot(plot_b, b2, '.', ms=15, c='k')
        axarr[1, col].annotate(r'$R^2$' + '=%0.4f' % r_value ** 2, xy=(0.02, 0.98), xycoords='axes fraction',
                               va='top', fontsize=20)

        axarr[1, col].tick_params(axis='both', labelsize=tick_size)
        if col != 0:
            axarr[1, col].get_xaxis().set_visible(False)
            axarr[1, col].get_yaxis().set_visible(False)
        else:
            axarr[1, col].locator_params(nbins=4)

    plt.tight_layout(h_pad=2, w_pad=2)
    plt.savefig('granger_figure.pdf', format='pdf')
