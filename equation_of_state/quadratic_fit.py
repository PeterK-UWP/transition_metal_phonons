import numpy as np


def quadratic_fit(array):
    x_values = array[0, :]  # first eneries[0] from two columns
    y_values = array[1, :]  # second eneries[1] from two columns

    quadratic_coefficients = np.polyfit(x_values, y_values, 2)

    return quadratic_coefficients