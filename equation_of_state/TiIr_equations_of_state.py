from two_column_text_read import two_column_text_read
from bivariate_statistics import bivariate_statistics
from quadratic_fit import quadratic_fit
from fit_curve_array import fit_curve_array
from equations_of_state import fit_eos
from convert_units import convert_units
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from format_plot import format_plot


import numpy as np
#\overline{3}
symmetries = {'Pm-3m': {'point_color': 'C0', 'data_point_shape': '.', 'line_style': 'solid', 'horizontal_offset': '', 'label': 'Pm3\u0305m'},
              'Im': {'point_color': 'k', 'data_point_shape': 'o', 'line_style': 'dotted', 'horizontal_offset': 'right', 'label': 'Im'},
              'P4mmm': {'point_color': 'g', 'data_point_shape': '.', 'line_style': 'dashed', 'horizontal_offset': 'left', 'label': 'P4/mmm'}}

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
             label=dictionary['label'])

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
plt.xlabel(r'$V$ [$\mathrm{\AA}^3$/f.u.]', fontsize=18)
plt.ylabel(r'$E$ [eV/f.u.]', fontsize=18)
plt.title('Equation of State TiIr', fontsize=20)
plt.legend()
plt.tight_layout()
plt.grid()
#plt.show()
plt.savefig("TiIrMurnaghanEquationsOfState.eps", bbox_inches="tight")

# cohesive energy, equilibrium volume, bulk modulus
