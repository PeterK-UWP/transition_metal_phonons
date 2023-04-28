import numpy as np
import re

"""
read file in
parse in the matrix from command prompt
calculate the new celldm()
output: celldm(1) = #
        celldm(2) = #
        celldm(3) = #
        celldm(5) = #
        
CELL_PARAMETERS (alat=  5.65028300)
   0.971011777   0.000000000   0.002016214
   0.000000000   0.970395586   0.000000000
  -0.002402789   0.000000000   1.159072501
"""


def calculate_celldm():
    lattice_parameter = 5.65028300
    parameters_matrix = np.array([
        [0.971011777, 0.000000000, 0.002016214],
        [0.000000000, 0.970395586, 0.000000000],
        [- 0.002402789, 0.000000000, 1.159072501]
    ])
    new_celldm_1 = lattice_parameter * parameters_matrix[0][0]
    new_celldm_2 = parameters_matrix[1][1]/parameters_matrix[0][0]
    new_celldm_5 = parameters_matrix[2][0]/np.sqrt(parameters_matrix[2][0]**2 + parameters_matrix[2][2]**2)
    new_celldm_3 = parameters_matrix[2][0]/(parameters_matrix[0][0]*new_celldm_5)
    return [new_celldm_1, new_celldm_2, new_celldm_3, new_celldm_5]

celldms = calculate_celldm()

for index, cell in enumerate(celldms):
    if index == 3:
        print(f'   celldm({index + 2}) = {cell:10.7f}')
    else:
        print(f'   celldm({index + 1}) = {cell:10.7f}')
