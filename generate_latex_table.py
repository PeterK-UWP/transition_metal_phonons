def generate_latex_table(quantities):
    number_of_keys = len(quantities)

    latex_table = "\\begin{tabular}{" + number_of_keys*'l' + "}\n"
    latex_table += f'{"Quantity":17}'

    for symmetry in quantities.keys():
        symmetry = rf"$" + f'{symmetry}' + rf"$"
        latex_table += rf"& {symmetry:>10} "
    latex_table += "\\\\\n"
    latex_table += "\\hline\n"
    first_symmetry = next(iter(quantities))
    for quantity in quantities[first_symmetry].keys():
        latex_table += f'{quantities[first_symmetry][quantity]["symbol"]:17}'
        for current_symmetry in quantities.keys():
            print(f' current_sym ={current_symmetry}, quan={quantity}, value={quantities[current_symmetry][quantity]["value"]:10}')
            latex_table += f'& {quantities[current_symmetry][quantity]["value"]:10} '
        latex_table += "\\\\\n"

    latex_table += "\\end{tabular}"
    return latex_table


if __name__ == "__main__":
    from equation_of_state.two_column_text_read import two_column_text_read
    from equation_of_state.convert_units import convert_units
    from equation_of_state.quadratic_fit import quadratic_fit
    from equation_of_state.equations_of_state import fit_eos
    from generate_dictionary import generate_dictionary
    from cohesive_energy import calculate_cohesive_energy

    import numpy as np

    # test_symmetries = ['Im', 'Pm-3m', 'P4mmm']
    # test_values = np.zeros(len(test_symmetries)*4)
    # test_dictionary = generate_dictionary(test_symmetries, test_values)
    # output_table = generate_latex_table(test_dictionary)
    # print()
    # print(output_table)

    symmetries = []
    values = []

    for symmetry in ['Im', 'Pm-3m', 'P4mmm']:
        filename = "TiIr." + symmetry + ".PBEsol.volume_totalenergy.txt"
        number_of_atoms = {"Ti": 1, "Ir": 1}
        # cut at 150
        energy_cutoff = {"Ti": -118.69491896*27.211396/2, "Ir": -181.04649751*27.211396/2}
        # check cutoff, ecut 140-Ti never converged, both should be the same cutoff

        directory = "data"

        volumes_energies = two_column_text_read(f'{directory}/{filename}')

        volumes_energies = np.array([convert_units(volumes_energies[0], 'cb/a'),
                                     convert_units(volumes_energies[1], 'rb/a')])

        quadratic_coefficients = quadratic_fit(volumes_energies)

        equation_of_state, equation_of_state_parameters = fit_eos(volumes_energies[0], volumes_energies[1],
                                                                  quadratic_coefficients)

        symmetries.append(symmetry)

        # Need to calculate E_coh from E_0 and order to V_0, E_coh, K_0, K_0'
        cohesive_energy = calculate_cohesive_energy(equation_of_state_parameters[0], energy_cutoff, number_of_atoms)
        print(cohesive_energy)
        for parameter in equation_of_state_parameters:
            values.append(parameter)
    print(symmetries, values)  # values array is correct
    eos_dictionary = generate_dictionary(symmetries, values)
    print(eos_dictionary)
    eos_latex_table = generate_latex_table(eos_dictionary)
    print(eos_latex_table)


# eq vol from each phase -> cubic angrtrom /fromula unit
# cohesive energy <- filesname to get total energy