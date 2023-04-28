from two_column_text_read import two_column_text_read
from bivariate_statistics import bivariate_statistics
from quadratic_fit import quadratic_fit
from fit_curve_array import fit_curve_array
# from plot_data_with_fit import plot_data_with_fit
from equations_of_state import fit_eos
from convert_units import convert_units
from numpy import linspace
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
display_graph = False


def parse_file_name(file_name):
    to_parse = file_name.split(".")
    symbol = to_parse[0]
    structure = to_parse[1]
    acronym = to_parse[2]
    return symbol, structure, acronym  # 1


file_name = "../data/Ti.P6_3mmc.PBEsol.volume_totalenergy.txt"
symbol, structure, acronym = parse_file_name(file_name)
array = two_column_text_read(file_name)  # 2
statistics = bivariate_statistics(array)  # 4
quadratic_coefficients = quadratic_fit(array)  # 5
# quadratic_coefficients = [quadratic_coefficients[0], quadratic_coefficients[1]]
print(quadratic_coefficients)

# print(array)
# mean_of_y, standard_deviation_of_y, x_min, x_max, y_min, y_max
# print(statistics)
min_x = statistics[2]
max_x = statistics[3]
"""
this comment provides the original graph I created using parameters before annotated graph. Un-comment and run to get
this graph. Save by commenting plt.show() out and uncomment #plt.savefig("Initial_plot.png") 
in plot_data_with_fit.py

fit_curve = fit_curve_array(quadratic_coefficients, min_x, max_x, number_of_points=100)
scatter_plot, curve_plot = plot_data_with_fit(array, fit_curve, data_format="bo", fit_format="k")
#plt.show()
plt.savefig("Initial_plot.png")



\/\/attempting the annotation \/\/
"""

undo_array = zip(*array)
array_2 = array  # list(undo_array)
# fit_eos_curve, bulk_modulus = fit_eos(array_2[0], array_2[1], quadratic_coefficients, eos='murnaghan', number_of_points=50)   #6
fit_eos_curve, fit_parameters = fit_eos(array[0], array[1], quadratic_coefficients, eos='murnaghan',
                                        number_of_points=50)  # 6
#plt.plot(linspace(array[0][0], array[0][-1], num=len(fit_eos_curve)), fit_eos_curve)
#plt.show()
print(fit_parameters)
bulk_modulus = fit_parameters[1]
equilibrium_volume = fit_parameters[3]


# print(fit_eos, bulk_modulus)

def annotate_graph(symbol, structure):
    ax.annotate(symbol, xy=(130, 0.001))

    ax.annotate(r'$ {}\overline{{{}}} {}$'.format(structure[0:2],
                                                  structure[3],
                                                  structure[1]),
                xy=(115, 0))

    ax.annotate('K_0={:.6f}GPa'.format(bulk_modulus_gpa),
                xy=(115, 0.001))

    ax.annotate('V_0={:.3f}A^3/atom'.format(eq_vol),
                xy=(115, -0.001))
    #plt.axvline(eq_vol - array_2[0][array_2[1].index(min(array_2[1]))] * 0.01, color="black", linestyle='--')

    plt.text(91, -0.0025, "created by Peter Kveton May/12/21")
    # plt.title("{} Equation of State for {} in DFT {}".format('Murnaghan', symbol, acronym))
    return ax, plt


fig = plt.figure()
ax = fig.add_subplot(111)
y_formatter = ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
volumes = linspace(min(array_2[0]), max(array_2[0]), len(fit_eos_curve))
plt.plot(array_2[0], array_2[1], 'o')
plt.plot(volumes, fit_eos_curve, color="black")

x_min = (min(array_2[0]) - (min(array_2[0]) * 0.10))
x_max = (max(array_2[0]) + (max(array_2[0]) * 0.10))
y_min = (-0.003)  # (min(array_2[1]) - (min(array_2[0]) * 0.00010))
y_max = (0.003)  # (max(array_2[1]) + (max(array_2[0]) * 0.00010))

plt.xlim(x_min, x_max)
#plt.ylim(y_min, y_max)
plt.xlabel(r'$V$ (Å$^3$/atom)')
plt.ylabel(r'$E$ (eV/atom)')
bulk_modulus_gpa = convert_units(bulk_modulus, "rb/cb")  # 7
eq_vol = fit_parameters[3]  # array_2[0][array_2[1].index(min(array_2[1]))]
array[1] = array[1] - fit_parameters[0]

#annotate_graph(symbol, structure)

fit_curve = fit_curve_array(quadratic_coefficients, min_x, max_x, number_of_points=100)
fit_curve[1] = fit_curve[1] - fit_parameters[0]
#print(array, fit_curve)
#scatter_plot, curve_plot = plot_data_with_fit(array, fit_curve, data_format="bo", fit_format="k")

plt.tight_layout()
if display_graph:
    plt.show()  # fix something above^^
elif not display_graph:
    plt.savefig("Ti.P6_3mmc.PBEsol.volume_total_energy.murnaghan_eos.png")
