import numpy as np


def fit_curve_array(quadratic_coefficients, min_x, max_x, number_of_points=100):

  if max_x < min_x:   #max values cannot be less than min values
    raise ArithmeticError
  if number_of_points <= 2:   #need more than 2 points
    raise IndexError

  x_array = np.linspace(min_x, max_x, number_of_points) #2 parameters
  y_array = np.polyval(quadratic_coefficients, x_array) #
  fit_curve = np.array((x_array, y_array))  #np.column_stack((x_array, y_array))

  return fit_curve