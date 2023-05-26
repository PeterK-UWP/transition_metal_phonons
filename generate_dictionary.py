def generate_dictionary(symmetries, values,
                        quantities=['equilibrium_volume', 'cohesive_energy', 'bulk_modulus', 'bulk_modulus_derivative'],
                        units=[r'$\AA^3$/f.u.', 'eV/f.u.', 'GPa', ' '],
                        symbols=[r'$V_0$', r'$E_\textrm{coh}$', r'$K_0$', r'$K^\prime_0$']):

    nested_dictionary_template = {}
    for quantity, unit, symbol in zip(quantities, units, symbols):
        nested_dictionary_template[quantity] = {'units': unit, 'symbol': symbol}

    dictionary = {}
    value_index = 0
    for symmetry in symmetries:
        dictionary[symmetry] = nested_dictionary_template
        for quantity in quantities:
            dictionary[symmetry][quantity]['value'] = values[value_index]
            value_index += 1
    return dictionary


if __name__ == '__main__':
    import numpy as np
    test_symmetries = ['Im', 'Pm3m', 'P4mmm', 'Cmmm']
    test_values = np.zeros(len(test_symmetries)*4)

    test_dictionary = generate_dictionary(test_symmetries, test_values)
    print(test_dictionary)
