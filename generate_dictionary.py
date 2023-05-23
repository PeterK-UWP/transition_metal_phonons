symmetries = ['Im', 'Pm3m', 'P4mmm', 'Cmmm']

quantities = ['equilibrium_volume', 'cohesive_energy', 'bulk_modulus', 'bulk_modulus_derivaitve']

units = ['$\AA^3$/f.u.', 'eV/f.u.', 'GPa', ' ']

symbols = ['V_0', 'E_\textrm{coh}', 'K_0', 'K^\prime_0']

dictionary = {}
for symmetry in symmetries:
    dictionary[symmetry] =