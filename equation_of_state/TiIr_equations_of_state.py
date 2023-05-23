from two_column_text_read import two_column_text_read
from bivariate_statistics import bivariate_statistics
from quadratic_fit import quadratic_fit
from fit_curve_array import fit_curve_array
from equations_of_state import fit_eos
from convert_units import convert_units
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

symmetries = {'Pm-3m': {'point_color': 'C0', 'data_point_shape': '.', 'line_style': 'solid', 'horizontal_offset': ''},
              'Im': {'point_color': 'C1', 'data_point_shape': '+', 'line_style': 'dotted', 'horizontal_offset': 'right'},
              'P4mmm': {'point_color': 'C2', 'data_point_shape': 'x', 'line_style': 'dashed', 'horizontal_offset': 'left'}}

for symmetry, dictionary in symmetries.items():
    filename = "TiIr." + symmetry + ".PBEsol.volume_totalenergy.txt"
    directory = "../data"

    volumes_energies = two_column_text_read(f'{directory}/{filename}')

    volumes_energies = np.array([convert_units(volumes_energies[0], 'cb/a'),
                                 convert_units(volumes_energies[1], 'rb/a')])

    statistics = bivariate_statistics(volumes_energies)

    quadratic_coefficients = quadratic_fit(volumes_energies)

    equation_of_state, equation_of_state_parameters = fit_eos(volumes_energies[0], volumes_energies[1],
                                                              quadratic_coefficients)
    print(equation_of_state_parameters)
    print(volumes_energies)
    volumes = np.linspace(volumes_energies[0][0], volumes_energies[0][-1], len(equation_of_state))
    plt.plot(volumes, equation_of_state, linestyle=dictionary['line_style'], color=dictionary['point_color'])
    plt.plot(volumes_energies[0], volumes_energies[1], dictionary['data_point_shape'], color=dictionary['point_color'],
             label=symmetry)
    """ 
    volume_range = volumes[-1] - volumes[0]
    match dictionary['horizontal_offset']:
        case 'right':
            plt.text(equation_of_state_parameters[-1] - 0.05 * volume_range, 0.5 * (volumes_energies[1][0] + volumes_energies[1][-1]), symmetry,
                     ha='right', va='center')
        case 'left':
            plt.text(equation_of_state_parameters[-1] + 0.05 * volume_range, 0.5 * (volumes_energies[1][0] + volumes_energies[1][-1]), symmetry,
                     ha='left', va='center')
        case _:
            plt.text(equation_of_state_parameters[-1], 0.5 * (volumes_energies[1][0] + volumes_energies[1][-1]), symmetry,
                     ha='center', va='center')
    """
plt.xlabel(r'$V$ [$\mathrm{\AA}^3$]')
plt.ylabel(r'$E$ [eV]')
plt.legend()
plt.tight_layout()
plt.show()

# cohesive energy, equilibrium volume, bulk modulus
