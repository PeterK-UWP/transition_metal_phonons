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


def parse_vc_relax_output_file(output_file):
    cell_parameters_lines = []
    number_of_lines_in_group = 4
    start_gather_string = 'CELL_PARAMETERS'
    stop_gather_string = 'ATOMIC_POSITIONS'
    with open(output_file, 'r') as file:
        gather = False
        for line in file:
            line = line.strip()
            if line.startswith(start_gather_string):
                gather = True
                cell_parameters_lines.append(line)
            elif line.startswith(stop_gather_string):
                gather = False
            elif line and gather:
                cell_parameters_lines.append(line)

    cell_parameters_lines = cell_parameters_lines[-number_of_lines_in_group:]
    lattice_constant = float(cell_parameters_lines[0].split()[2][:-1])
    lattice_vectors = np.empty([3, 3])
    for index, line in enumerate(cell_parameters_lines[1:]):
        lattice_vectors[index] = [float(entry) for entry in line.split()]
    return lattice_constant, lattice_vectors


def calculate_celldm(lattice_parameter, parameters_matrix):
    # lattice_parameter = 5.65028300
    # parameters_matrix = np.array([
    #     [0.971011777, 0.000000000, 0.002016214],
    #     [0.000000000, 0.970395586, 0.000000000],
    #     [- 0.002402789, 0.000000000, 1.159072501]
    # ])
    new_celldm_1 = lattice_parameter * parameters_matrix[0][0]
    new_celldm_2 = parameters_matrix[1][1]/parameters_matrix[0][0]
    new_celldm_5 = parameters_matrix[2][0]/np.sqrt(parameters_matrix[2][0]**2 + parameters_matrix[2][2]**2)
    new_celldm_3 = parameters_matrix[2][0]/(parameters_matrix[0][0]*new_celldm_5)
    return [new_celldm_1, new_celldm_2, new_celldm_3, new_celldm_5]


def main():
    vc_relax_output_file = 'data/TiIr.Im.vc-relax.out'
    lattice_parameter, lattice_matrix = parse_vc_relax_output_file(vc_relax_output_file)
    celldms = calculate_celldm(lattice_parameter, lattice_matrix)

    for index, cell in enumerate(celldms):
        if index == 3:
            print(f'   celldm({index + 2}) = {cell:10.7f}')
        else:
            print(f'   celldm({index + 1}) = {cell:10.7f}')
    return


if __name__ == '__main__':
    main()

