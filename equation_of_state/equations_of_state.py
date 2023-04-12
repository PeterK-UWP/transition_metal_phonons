import numpy as np


def fit_eos(volumes, energies, quadratic_coefficients, eos='vinet', number_of_points=50):
    """
    Returns a NumPy array of values evaluated to an equation of state fit
    :param volumes:                 NumPy array(N) :: volumes (x-values) to be fit
    :param energies:                NumPy array(N) :: energies (y-values) to be fit
    :param quadratic_coefficients:  list(3) :: coefficients of the quadratic polynomial already fit to the data
    :param eos:                     str :: equation of state name ('murnaghan', 'birch-murnaghan', 'vinet')
    :param number_of_points:        int, optional :: number of points to evaluate fit function on
    :return:                        NumPy array(number_of_points) :: equation of state fit evaluated on grid,

    starting at volumes[0] and ending at volume[-1]
    """
    from scipy.optimize import curve_fit

    # Dictionary holding lambda functions from current module.
    lambda_dictionary = {
        'vinet': globals()['vinet'],
        'murnaghan': globals()['murnaghan'],
        'birch-murnaghan': globals()['birch_murnaghan']
    }

    # Get extremes of data and calculate range

    minimum_volume = np.amin(volumes)
    maximum_volume = np.amax(volumes)

    #   axis of symmetry: x = -b / 2a
    quadratic_axis_of_symmetry = -quadratic_coefficients[1] / (2 * quadratic_coefficients[0])
    #   minimum: y = -b^2 / (4 a)  + c
    quadratic_minimum = -quadratic_coefficients[1] ** 2 / (4 * quadratic_coefficients[0]) + quadratic_coefficients[2]
    #   bulk modulus: K_0 = 2 * a / V_0 for E(V) = a*V^2 + b*V + E0
    quadratic_bulk_modulus = 2. * quadratic_coefficients[0] / quadratic_axis_of_symmetry

    bulk_modulus_derivative = 3.7

    # Get realistic equation of state fit

    initial_parameters = [quadratic_minimum, quadratic_bulk_modulus,
                          bulk_modulus_derivative, quadratic_axis_of_symmetry]

    eos_parameters, eos_covariances = curve_fit(lambda_dictionary[eos.lower()], volumes, energies,
                                                p0=initial_parameters)
    fit_curve_volumes = np.linspace(minimum_volume, maximum_volume, num=number_of_points)
    eos_fit_curve = lambda_dictionary[eos.lower()](fit_curve_volumes,
                                                   eos_parameters[0], eos_parameters[1], eos_parameters[2],
                                                   eos_parameters[3])

    return eos_fit_curve, eos_parameters


def murnaghan(volumes, equilibrium_energy, bulk_modulus, bulk_modulus_derivative, equilibrium_volume):
    """
    Murnaghan equation of state: E(V) = E_0 + K_0 V_0 [ (1 / (K_0' (K_0' - 1))) (V / V_0)^(-(K_0' - 1)) +
                                                        (1 / K_0') (V / V_0) -

                                                        (1 / (K_0' - 1)) ]

    :param volumes:                 NumPy array of volumes per atom

    :param equilibrium_energy:      equilibrium energy E_0

    :param bulk_modulus:            bulk modulus K_0 = ∂^E/(∂V)^2 / V_0

    :param bulk_modulus_derivative: pressure derivative of bulk modulus ∂K_0/∂P

    :param equilibrium_volume:      equilibrium volume V_0

    :return:                        NumPy array of Murnaghan equation of state values at input volumes
    """
    k0pm1 = bulk_modulus_derivative - 1.0  # K_0' - 1
    return equilibrium_energy + (bulk_modulus * equilibrium_volume *
                                 (((1.0 / (bulk_modulus_derivative * k0pm1)) *
                                   np.power((volumes / equilibrium_volume), (-k0pm1))) +
                                  (volumes / (bulk_modulus_derivative * equilibrium_volume)) - (1.0 / k0pm1)))


def birch_murnaghan(volumes, equilibrium_energy, bulk_modulus, bulk_modulus_derivative, equilibrium_volume):
    """
    Birch-Murnaghan equation of state: E(V) = E_0 + (9/16) K_0 V_0 {[ (V / V_0)^(-(2/3)) - 1 ]^3 K_0' +
                                                                    [ (V / V_0)^(-(2/3)) - 1]^2 *
                                                                    [ 6 - 4 (V / V_0)^(-(2/3)) ]}

    :param volumes:                 NumPy array of volumes per atom

    :param equilibrium_energy:      equilibrium energy E_0

    :param bulk_modulus:            bulk modulus K_0 = ∂^E/(∂V)^2 / V_0

    :param bulk_modulus_derivative: pressure derivative of bulk modulus ∂K_0/∂P

    :param equilibrium_volume:      equilibrium volume V_0

    :return:                        NumPy array of the Birch-Murnaghan equation of state values at input volumes
    """
    reduced_volume_area = np.power(volumes / equilibrium_volume, -2. / 3.)
    return equilibrium_energy + (9. * bulk_modulus * equilibrium_volume / 16.) * (
            np.power(reduced_volume_area - 1., 3.) * bulk_modulus_derivative +
            np.power(reduced_volume_area - 1., 2.) * (6. - 4. * reduced_volume_area))


def vinet(volumes, equilibrium_energy, bulk_modulus, bulk_modulus_derivative, equilibrium_volume):
    """
    Vinet equation of state: E(V) = E_0 + (2 K_0 V_0 / (K_0' - 1)^2) *
                                        - {2 - [5 + 3 (V / V_0)^(1/3) (K_0' - 1) - 3 K_0']
                                               exp(- (3/2) (K_0' - 1) (1 - (V / V_0)^(1/3))})

    :param volumes:                 NumPy array of volumes per atom

    :param equilibrium_energy:      equilibrium energy E_0

    :param bulk_modulus:            bulk modulus K_0 = ∂^E/(∂V)^2 / V_0

    :param bulk_modulus_derivative: pressure derivative of bulk modulus ∂K_0/∂P

    :param equilibrium_volume:      equilibrium volume V_0

    :return: NumPy array of the Vinet equation of state values at input volumes
    """
    k0pm1 = bulk_modulus_derivative - 1  # K_0' - 1
    k0pm1_squared = np.power(k0pm1, 2)
    reduced_volume_lengths = np.cbrt(volumes / equilibrium_volume)

    exponential_argument = -1.5 * k0pm1 * (reduced_volume_lengths - 1.)
    try:
        exponential_factor = np.exp(exponential_argument)
        vinet_eos = equilibrium_energy + \
                    (2. * bulk_modulus * equilibrium_volume / k0pm1_squared) * \
                    (2. - (5. + 3. * reduced_volume_lengths * k0pm1 - 3 * bulk_modulus_derivative) *
                     exponential_factor)
    except:
        # assuming the failure is factur
        vinet_eos = equilibrium_energy + \
                    (2. * bulk_modulus * equilibrium_volume / k0pm1_squared) * \
                    (2. - (5. + 3. * reduced_volume_lengths * k0pm1 - 3 * bulk_modulus_derivative))

    return vinet_eos
