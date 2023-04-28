from matplotlib import rcParams
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants


# THz_per_inverse_cm = 0.02998
# (10^-12 THz/Hz) (10^-2 cm / m)
THz_per_inverse_cm = constants.value(u'inverse meter-hertz relationship')*(10**-12)/(10**-2)
band_color = 'blue'


def parse_command_line_arguments():
    """
    Take matdyn input file name from command line if provided
    :return: matdyn input filename as string
    """
    import argparse
    parser = argparse.ArgumentParser(description='Process matdyn.x phonon frequencies...')
    parser.add_argument('-i', metavar="(matdyn.in)", type=str,
                        help='name of the matdyn.x input file', default='matdyn.in')
    parser.add_argument('--verbose', '-v', action=argparse.BooleanOptionalAction)
    parser.add_argument('--save_plot', '-s', action=argparse.BooleanOptionalAction)
    arguments = parser.parse_args()
    return arguments.i, arguments.verbose, arguments.save_plot


def get_frequencies_filename(input_filename):
    """
    Look for name of frequencies file in the matdyn input file
    :param input_filename:
    :return: frequencies filename as string
    """
    filename = 'matdyn.freq'
    with open(input_filename) as matdyn_file:
        for line in matdyn_file:
            line = line.strip()
            if line[:5] == 'flfrq':
                filename = line[7:len(line)-1]

    return filename


def extract_input_file_qpoints(filename):
    points = []
    point_numbers = []
    point_labels = []
    with open(filename) as matdyn_file:
        for line in matdyn_file:
            if line.count("%") > 0:
                line = line.split()
                points.append([float(line[0]), float(line[1]), float(line[2])])
                point_numbers.append(int(line[-2]))
                if line[-1][1:] == "Gamma":
                    point_labels.append(r'$\Gamma$')
                else:
                    point_labels.append(line[-1][1:])

    return points, point_numbers, point_labels


def calculate_high_symmetry_q_point_indices(point_numbers):
    point_indices = [0]
    index = 0
    for number in point_numbers[:-1]:
        index += number
        point_indices.append(index)

    return point_indices


def parse_frequencies_file(filename):
    import re
    with open(filename) as frequency_file:
        first_line = frequency_file.readline()
        number_of_bands = int(re.sub('\,', '', first_line).split()[2])
        number_of_points = int(first_line.split()[4])
        # print(f"{number_of_bands} bands with {number_of_points} q-points found in {filename}")
        number_of_frequencies_per_line = 6  # matdyn.x outputs six values per line
        if number_of_bands >= number_of_frequencies_per_line:
            point_line_modulus = int(number_of_bands / number_of_frequencies_per_line) + 1
        else:
            point_line_modulus = 2
        # print(f"Expecting q-point value on lines of index {point_line_modulus}")
        frequency_values = np.zeros((number_of_bands, number_of_points))

        frequency_line_index = 0
        for line_index, line in enumerate(frequency_file.readlines()):
            if line_index % point_line_modulus != 0:
                point_index = int(line_index/point_line_modulus)
                frequency_values[frequency_line_index:frequency_line_index + number_of_frequencies_per_line,
                                 point_index] = line.split()
                # print(line.split())
                frequency_line_index += number_of_frequencies_per_line
            else:
                frequency_line_index = 0

    return frequency_values


def get_q_point_plot_values(number_of_points, filename):
    point_plot_values = np.zeros(number_of_points)
    with open(filename) as frequencies_gnuplot_file:
        for line_index, line in enumerate(frequencies_gnuplot_file.readlines()):
            point_plot_values[line_index] = line.split()[0]

    return point_plot_values


def get_high_symmetry_q_point_plot_values(point_plot_values, point_indices):
    high_symmetry_plot_values = []
    for point_index in point_indices:
        high_symmetry_plot_values.append(point_plot_values[point_index])
    return high_symmetry_plot_values


def plot_frequency_bands(high_symmetry_plot_values, high_symmetry_labels, plot_values, frequency_bands,
                           bands_color='black', filename=''):

    # figure, axes = plt.subplots()
    rcParams['text.usetex'] = False

    if rcParams['text.usetex']:
        rcParams["axes.formatter.use_mathtext"] = True
        rcParams['font.sans-serif'] = "cmr12"

    # print(f'ω_min, ω_max = {np.amin(frequency_bands), np.amax(frequency_bands)} cm^-1')

    # Convert frequency units to THz
    frequency_bands *= THz_per_inverse_cm

    # Format horizontal axis
    plt.xlim([high_symmetry_plot_values[0], high_symmetry_plot_values[-1]])
    plt.xlabel(r'$\vec{q}$')
    plt.xticks(high_symmetry_plot_values, high_symmetry_labels)

    # Plot zero-frequency line
    plt.axhline(0., color='black')

    # Format vertical axis
    plt.ylim([np.floor(np.amin(frequency_bands)), np.ceil(np.amax(frequency_bands))])
    #plt.ylabel(r'$\omega(\vec{q})\textrm{ [THz]}$')


    # Mark high symmetry points vertically
    for plot_value in high_symmetry_plot_values:
        plt.axvline(plot_value, color='black', linestyle='--')

    # Make the bands plot
    number_of_bands = len(frequency_bands)
    for band_index in range(0, number_of_bands-1):
        plt.plot(plot_values, frequency_bands[band_index], color=bands_color)

    if filename == '':
        plt.show()
    else:
        figure_filename = filename[:-4] + 'png'
        plt.savefig(figure_filename)
    return


def parse_frequencies(input_parameters):
    if len(input_parameters) == 0:
        input_file, verbose, save_plot = parse_command_line_arguments()
    else:
        input_file = input_parameters[0]
        verbose = input_parameters[1]
        save_plot = input_parameters[2]

    if verbose:
        print(f"Matdyn input file to parse: {input_file}")

    frequencies_filename = get_frequencies_filename(input_file)
    if verbose:
        print(f"Frequencies file to parse: {frequencies_filename}")

    gnuplot_frequencies_filename = frequencies_filename + ".gp"
    if verbose:
        print(f"Frequencies gnuplot format: {gnuplot_frequencies_filename}")

    high_symmetry_q_points, q_point_numbers, high_symmetry_q_point_labels = extract_input_file_qpoints(input_file)
    number_of_high_symmetry_q_points = len(high_symmetry_q_point_labels)
    number_of_q_points = sum(q_point_numbers)
    if verbose:
        print(f"{number_of_high_symmetry_q_points} high-symmetry q-points found")
        print(f"\t{number_of_q_points} total q-points to be read in")

    high_symmetry_q_point_indices = calculate_high_symmetry_q_point_indices(q_point_numbers)
    if verbose:
        print(f"\t{high_symmetry_q_point_indices} indices found for high-symmetry q-points")

    frequencies = parse_frequencies_file(frequencies_filename)
    if verbose:
        print(f'ω_min, ω_max = {np.amin(frequencies), np.amax(frequencies)} cm^-1')
    number_of_bands = len(frequencies)
    if verbose:
        print(f"{number_of_bands*len(frequencies[0])} frequencies read in {number_of_bands} bands over "
              f"{len(frequencies[0])} q-points")

    q_point_plot_values = get_q_point_plot_values(number_of_q_points, gnuplot_frequencies_filename)
    high_symmetry_q_point_plot_values = get_high_symmetry_q_point_plot_values(q_point_plot_values,
                                                                          high_symmetry_q_point_indices)
    if verbose:
        for plot_value, point_label in zip(high_symmetry_q_point_plot_values, high_symmetry_q_point_labels):
            print(f"x-coordinate {plot_value} found for point {point_label}")

    if save_plot:
        plot_frequency_bands(high_symmetry_q_point_plot_values, high_symmetry_q_point_labels, q_point_plot_values,
                             frequencies, bands_color=band_color, filename=frequencies_filename)
    else:
        plot_frequency_bands(high_symmetry_q_point_plot_values, high_symmetry_q_point_labels, q_point_plot_values,
                             frequencies, bands_color=band_color)

    return


if __name__ == "__main__":
    import os
    os.chdir('data')

    matdyn_input_file = 'Ir.Fm-3m.matdyn.in'
    verbose_output = False
    make_plot_file = False

    parse_frequencies([matdyn_input_file, verbose_output, make_plot_file])



