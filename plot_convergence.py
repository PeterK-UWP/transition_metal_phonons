"""
determine file (string) GREAT!!
read data from file :D
convert data units (data in Rydbergs convert to ) we did it!!
determine convergence threshold :)
plot data Congratulations!!
mark/label convergence threshold
"""
import numpy as np
import matplotlib.pyplot as plt
from conversions import rydberg_to_eV

file_name = 'Ir.Fm-3m.PBEsol.ecutwfc_totalenergy.txt'
convergence_threshold = 0.001  # eV
data = np.loadtxt(file_name)
data_labels = file_name.split(".")[3].split("_")
# print(f'{data_labels}')
# print(data.transpose())
x_values = data.transpose()[0]
y_values = data.transpose()[1]

if data_labels[1] == 'totalenergy':
    y_values = y_values * rydberg_to_eV

y_value_differences = y_values - y_values[-1]
# np.abs(np.diff(y_values))  # apply later
y_value_threshold_differences = np.abs(np.abs(y_value_differences) - convergence_threshold)
converged_index = np.argmin(y_value_threshold_differences) + 1

if (np.abs(y_value_differences) - convergence_threshold)[converged_index] < 0:   # error
    converged_index -= 1

converged_y_value = y_values[converged_index]

plt.plot(x_values, y_values)
plt.axhline(converged_y_value, linestyle="-")
if data_labels[0] == 'ecutwfc':
    x_label = r'$E_\mathrm{cut}$ (Ry)'
if data_labels[1] == 'totalenergy':
    y_label = r'$E_\mathrm{total}$ (eV)'
y_minimum = np.min(y_values) - convergence_threshold
y_maximum = y_minimum + (10 * convergence_threshold)
# plt.ylim([y_minimum, y_maximum])
plt.xlim([np.min(x_values[-9]), x_values[-1]])
plt.ylim([np.min(y_values), y_values[-9]])
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.show()
"""
Need more cutoff energies for Ir Fm-3m  :(
"""
