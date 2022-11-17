def plot_convergence(file_name, convergence_threshold, fit_start_index,
                     write_figure_to_file=True, display_figure=True):
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
    from scipy.optimize import curve_fit
    from pathlib import Path

    data = np.loadtxt(file_name)
    data_labels = file_name.split(".")[3].split("_")
    # print(f'{data_labels}')
    # print(data.transpose())
    x_values = data.transpose()[0]
    y_values = data.transpose()[1]

    atomic_file_name = file_name.split(".")[0] + '.Pm-3m'
    for name_component in file_name.split(".")[2:]:
        atomic_file_name += '.' + name_component
    # print(atomic_file_name)

    path = Path(atomic_file_name)
    if path.is_file():
        atomic_data = np.loadtxt(atomic_file_name)
        atomic_x_values = atomic_data.transpose()[0]
        atomic_y_values = atomic_data.transpose()[1]
        Delta_x_values = []
        Delta_y_values = []
        for index, x_value in enumerate(x_values):
            atomic_index = np.where(atomic_x_values == x_value)[0]
            if len(atomic_index) > 0:
                Delta_x_values.append(x_value)
                Delta_y_values.append(y_values[index] - atomic_y_values[atomic_index[0]])
        x_values = np.array(Delta_x_values)
        y_values = np.array(Delta_y_values)

    # Convert units
    if data_labels[1] == 'totalenergy':
        y_values = y_values * rydberg_to_eV

    # Identify first y-value within convergence threshold
    # Δy_i = y_i - y_N
    y_value_differences = y_values - y_values[-1]
    # |Δy|_i
    y_value_absolute_differences = np.abs(y_value_differences)
    # | |Δy|_i - Δy_convergence |
    y_value_threshold_differences = np.abs(y_value_absolute_differences - convergence_threshold)
    # i_converged = i_min(|Δy| - Δy_convergence) + 1
    converged_index = np.argmin(y_value_threshold_differences) + 1

    # In the case of no convergence, set the converged index to the last entry in the y_values array
    if converged_index == len(y_values):
        converged_index -= 1

    # if (np.abs(y_value_differences) - convergence_threshold)[converged_index] < 0:   # error
    #    converged_index -= 1

    # y_i_converged
    converged_y_value = y_values[converged_index]

    # Determine axis limits
    #   Set y values up by convergence threshold intervals above and below the minimum value
    y_minimum = (y_values[-1]) - 2 * convergence_threshold
    y_maximum = (y_values[-1]) + (8 * convergence_threshold)
    plt.ylim([y_minimum, y_maximum])
    #   Set x values by lowest y-value appearing off the chart
    off_chart_index = 0
    for y_value in y_values:
        if y_value > y_maximum:
            off_chart_index += 1
    off_chart_index -= 1
    x_minimum = x_values[off_chart_index]
    x_range = x_values[-1] - x_minimum
    x_maximum = x_minimum + 1.1 * x_range
    # plt.xlim([x_minimum, x_maximum])

    # Plot values
    plt.plot(x_values, y_values, 'o')

    # Fit values
    def fit_function(x, c0, c1, c2, c3):
        # y = d + a e^(b (x-c)) / (x-c) where b is negative
        # return c0 * np.exp(c1*(x - c2)) * np.power(x - c2, -1) + c3
        #
        # y = d + a b / (x - c) where b is redundant
        # return c0 * c1 * np.power(x - c2, -1) + c3
        #
        # y = d + a e^(-b (x-c)) / (x-c) where b is positive
        return c0 * np.exp(-c1 * (x - c2)) * np.power(x - c2, -1) + c3

    # Parameter guesses that work
    c3_guess = y_values[-1]
    c2_guess = 0
    c0_guess = 1
    c1_guess = 0

    # Parameter guesses that don't work
    # c2_guess = np.min(x_values[fit_start_index:])
    # logarithmic_argument = y_values[fit_start_index] - c3_guess
    # logarithmic_argument *= (x_values[fit_start_index] - c2_guess)/c0_guess
    # c1_guess = np.log(logarithmic_argument)/(x_values[fit_start_index] - c2_guess)

    parameter_guess = [c0_guess, c1_guess, c2_guess, c3_guess]
    print(f'parameters give to fit = {parameter_guess}')
    fit_parameters, fit_covariance = curve_fit(fit_function, x_values[fit_start_index:], y_values[fit_start_index:],
                                               p0=parameter_guess)
    print(f'parameters produced by fit = {fit_parameters}')
    fit_x_values = np.linspace(x_minimum, x_maximum)
    fit_y_values = fit_function(fit_x_values, *fit_parameters)
    plt.plot(fit_x_values, fit_y_values, color='orange')

    # Plot convergence threshold lines
    # plt.axhline(np.min(y_values) - convergence_threshold, linestyle="--", color='black')
    # plt.axhline(np.min(y_values) + convergence_threshold, linestyle="--", color='black')
    plt.axhline((y_values[-1]) - convergence_threshold, linestyle="--", color='black')
    plt.axhline((y_values[-1]) + convergence_threshold, linestyle="--", color='black')
    # Plot convergence point vertical line segment
    plt.axvline(x_values[converged_index],
                ymax=(y_values[converged_index] - y_minimum) / (y_maximum - y_minimum), linestyle="-.", color="black")

    # Label convergence thresholds
    upper_threshold_label = f'+{convergence_threshold:.3f} eV'
    lower_threshold_label = f'-{convergence_threshold:.3f} eV'

    plot_x_minimum, plot_x_maximum = plt.xlim()
    plot_x_range = plot_x_maximum - plot_x_minimum
    threshold_label_x = 0.03 * plot_x_range + plot_x_minimum

    plot_y_minimum, plot_y_maximum = plt.ylim()
    plot_y_range = plot_y_maximum - plot_y_minimum
    upper_threshold_label_y = 0.05 * plot_y_range + ((y_values[-1]) + convergence_threshold)
    lower_threshold_label_y = -0.05 * plot_y_range + ((y_values[-1]) - convergence_threshold)
    plt.text(threshold_label_x, upper_threshold_label_y, upper_threshold_label, va='center')
    plt.text(threshold_label_x, lower_threshold_label_y, lower_threshold_label, va='center')

    # Label axes
    if data_labels[0] == 'ecutwfc':
        x_label = r'$E_\mathrm{cut}$ (Ry)'
    if data_labels[1] == 'totalenergy':
        y_label = r'$E_\mathrm{total}$ (eV)'
        if path.is_file():
            y_label = r'$\Delta E_\mathrm{total}$ (eV)'

    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Label plot
    # chemical_label = file_name.split('.')[0].split('/')[1]
    chemical_label = file_name.split('.')[0]

    structure_label = file_name.split('.')[1]
    #    Convert "-" to overline
    if len(structure_label.split('-')) > 1:
        structure_label = r'$' + structure_label.split('-')[0] + '\overline{' + structure_label.split('-')[1][0] + '}' + \
                          structure_label.split('-')[1][1:] + '$'
    exchange_correlation_label = file_name.split('.')[2]
    plot_label = chemical_label + ' ' + structure_label + ' ' + exchange_correlation_label
    plot_label_x = plot_x_maximum - 0.05 * plot_x_range
    plot_label_y = plot_y_maximum - 0.10 * plot_y_range
    plt.text(plot_label_x, plot_label_y, plot_label, fontsize=16, ha='right')

    if write_figure_to_file:
        structure_label = file_name.split('.')[1]
        figure_file_name = chemical_label + '.' + structure_label + '.' + exchange_correlation_label + '.' + \
                           data_labels[0] + '_' + data_labels[1] + '.png'
        print(f'Writing figure to {figure_file_name}')
        plt.savefig(figure_file_name, format='png')

    if display_figure:
        plt.show()
    """
    Need more cutoff energies for Ir Fm-3m  :(
    """
    return


if __name__ == '__main__':
    # file_of_interest = 'data/Ir.Fm-3m.PBEsol.ecutwfc_totalenergy.txt'
    file_of_interest = 'Pt.Fm-3m.PBEsol.ecutwfc_totalenergy.txt'
    energy_convergence_threshold = 0.001  # eV
    start_index_for_fit = 15
    plot_convergence(file_of_interest, energy_convergence_threshold, start_index_for_fit,
                     write_figure_to_file=True, display_figure=True)