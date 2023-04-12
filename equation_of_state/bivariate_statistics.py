import numpy as np
from scipy import stats


def bivariate_statistics(data):
    if len(data) != 2 or len(data[0]) <= 1:
        raise IndexError

    stat = stats.stats.describe(data, axis=1)
    mean_of_y = stat.minmax[0][1]
    x_min, x_max = stat.minmax[0][0], stat.minmax[1][0]
    y_min, y_max = stat.minmax[0][1], stat.minmax[1][1]
    standard_deviation_of_y = np.sqrt(stat.variance[1])
    statistics = np.array([mean_of_y, standard_deviation_of_y, x_min, x_max, y_min, y_max])

    return statistics