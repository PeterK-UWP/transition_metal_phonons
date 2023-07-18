

def generate_dictionary(symmetries, values,
                        quantities=['cohesive_energy', 'bulk_modulus', 'bulk_modulus_derivative', 'equilibrium_volume'],
                        units=['eV/f.u.', 'GPa', ' ', r'$\AA^3$/f.u.'],
                        symbols=[r'$E_\textrm{coh}$', r'$K_0$', r'$K^\prime_0$', r'$V_0$']):
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

"""
def atom_energies_dictionary(symmetries):
    For TiIr eos:
        Im    Minimum ->  Volume: 197.0128      Energy: -300.93892729
        Pm3m  Minimum ->  Volume: 197.94760000  Energy: -300.92557383
        P4mmm Minimum ->  Volume: 196.9374      Energy: -300.93892925
    For Ti
        P63mmm: ecut: 150  energy: -238.29124857
        Pm3m:   ecut: 150  energy: -118.69722486
    For Ir
        Fm3m:   ecut: 140  energy: -181.66967817
        Pm3m:   ecut: ???  energy: ???     the plot didnt work... currently loaded
    return dictionary

def atom_numbers_dictionary(symmytries, molecular_energy):
    Ti: 1
    Ir: 1
    return dictionary
"""

if __name__ == '__main__':
    import numpy as np

    test_symmetries = ['Im', 'Pm3m', 'P4mmm']  # , 'Cmmm']

    test_values = np.arange(4 * len(test_symmetries))
    # test_values = np.zeros(len(test_symmetries)*4)
    print(test_values)

    test_dictionary = generate_dictionary(test_symmetries, test_values)
    print(test_dictionary)
