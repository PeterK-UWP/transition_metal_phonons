import numpy as np

"""
dict of isolated energies
dict of weights "number of atoms/formula unit"
"""


def calculate_cohesive_energy(condensed_energy, energy_dictionary, atoms):
    # isolated_energy = [number*energy for key, number in atoms.items() for energy in energy_dictionary.values()
    #                       if key in energy_dictionary]
    products = {atom: energy_dictionary[atom] * atoms[atom] for atom in energy_dictionary.keys() & atoms.keys()}
    isolated_energy = sum(products.values())
    cohesive_energy = condensed_energy - isolated_energy

    return cohesive_energy


if __name__ == "__main__":
    molecular_energy = -20
    atomic_energies = {"C": -10, "O": -4}
    number_of_atoms = {"C": 1, "O": 2}
    # molecular_energy = -118.69722486 + -181.04649751
    # atomic_energies = {"Ti": -118.69722486, "Ir": -181.04649751}  # ecut 150
    # number_of_atoms = {"Ti": 1, "Ir": 1}
    co2_cohesive_energy = calculate_cohesive_energy(molecular_energy, atomic_energies, number_of_atoms)
    print(co2_cohesive_energy)
