

def generate_dictionary(symmetries, values,
                        quantities=['equilibrium_volume', 'cohesive_energy', 'bulk_modulus', 'bulk_modulus_derivative'],
                        units=[r'$\AA^3$/f.u.', 'eV/f.u.', 'GPa', ' '],
                        symbols=[r'$V_0$', r'$E_\textrm{coh}$', r'$K_0$', r'$K^\prime_0$']):
    dictionary = {}
    value_index = 0

    for symmetry in symmetries:
        selected_values = values[value_index * len(quantities): value_index * len(quantities) + len(quantities)]
        symmetry_dictionary = {}
        for quantity, symbol, value, unit in zip(quantities, symbols, selected_values, units):
            symmetry_dictionary[quantity] = {'symbol': symbol,
                                             'value': value,
                                             'unit': unit}
        dictionary[symmetry] = symmetry_dictionary
        # for symbol, value, unit in zip(symbols, selected_values, units):
        value_index += 1
    return dictionary


if __name__ == '__main__':
    import numpy as np

    test_symmetries = ['Im', 'Pm3m', 'P4mmm']  # , 'Cmmm']

    test_values = np.arange(4 * len(test_symmetries))
    # test_values = np.zeros(len(test_symmetries)*4)
    print(test_values)

    test_dictionary = generate_dictionary(test_symmetries, test_values)
    print(test_dictionary)
