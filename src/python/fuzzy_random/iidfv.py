
from itertools import izip
import random

from scipy import comb
import numpy as np

from common import random_from_distribution, prepare_probabilities

__all__ = ["iidfv", "iidfv_discrete"]

def iidfv_fast(probabilities, number_of_variables):
    if number_of_variables <= 0:
        raise ValueError("NUMBER_OF_VARIABLES must be a positive number")

    if len(probabilities) <= 0:
        raise ValueError("PROBABILITIES must be nonempty")

    probabilities = np.array(sorted(probabilities, reverse=True))
    number_of_probabilities = len(probabilities)

    variables = np.zeros(shape=(number_of_variables,), dtype=int)

    coeffs = np.float64(np.vander(np.arange(number_of_probabilities), number_of_variables - 1))
    binom_coeffs = np.column_stack(
        comb(number_of_variables - k, np.arange(1, number_of_variables))
        for k in xrange(1, number_of_variables))
    binom_coeffs = binom_coeffs.T

    raw_counts = np.zeros(shape=(number_of_probabilities, number_of_variables - 1))
    for event in xrange(number_of_probabilities):
        for variable in xrange(number_of_variables - 1):
            raw_counts[event][variable] = np.dot(coeffs[event][variable:],
                                                 binom_coeffs[variable][:number_of_variables - 1 - variable])
    raw_counts = raw_counts.T

    distribution = np.zeros(shape=(number_of_probabilities,))
    for current_variable in xrange(number_of_variables):

        for current_event in xrange(number_of_probabilities):

            least_probable_event = max(variables[:current_variable].tolist() + [current_event])

            if current_variable < number_of_variables - 1:
                current_event_counts = np.zeros(shape=(number_of_probabilities,))
                current_event_counts[:least_probable_event] = 0
                current_event_counts[least_probable_event] = np.sum(raw_counts[current_variable][:least_probable_event + 1])
                current_event_counts[least_probable_event + 1:] = raw_counts[current_variable][least_probable_event + 1:]

                distribution[current_event] = np.dot(current_event_counts, probabilities)

            else:
                distribution[current_event] = probabilities[least_probable_event]

        distribution /= sum(distribution)
        variables[current_variable] = random_from_distribution(distribution)

    return variables

def iidfv_naive(probabilities, number_of_variables):
    if number_of_variables <= 0:
        raise ValueError("NUMBER_OF_VARIABLES must be a positive number")

    if len(probabilities) <= 0:
        raise ValueError("PROBABILITIES must be nonempty")

    probabilities = np.array(sorted(probabilities, reverse=True))
    number_of_probabilities = len(probabilities)

    variables = np.zeros(shape=(number_of_variables,))
    stack = np.zeros(shape=(number_of_variables,), dtype=int)
    distribution = np.zeros(shape=(number_of_probabilities,))

    for current_variable in xrange(number_of_variables):
        if current_variable > 0:
            stack[current_variable - 1] = variables[current_variable - 1]

        for current_event in xrange(number_of_probabilities):
            distribution[current_event] = 0

            stack[current_variable] = current_event
            min_beginning = min(probabilities[stack[:current_variable + 1]])
            stack_remainder_len = number_of_variables - (current_variable + 1)
            if stack_remainder_len > 0:
                for stack_remainder in product(xrange(number_of_probabilities),
                                               repeat=stack_remainder_len):
                    stack_remainder = np.array(stack_remainder)
                    min_remainder = min(probabilities[stack_remainder])
                    distribution[current_event] += min(min_beginning, min_remainder)
            else:
                distribution[current_event] = min_beginning

        distribution /= sum(distribution)
        variables[current_variable] = random_from_distribution(distribution)

    return variables

iidfv_discrete = iidfv_fast

def continuous_variable(variable, granules):
    granule = granules[variable]
    interval = random.choice(granule.intervals)
    value = random.uniform(*interval)
    if value in [float("+Inf"), float("-Inf")]:
        print value,
    return value

def iidfv(coordinates):
    if len(coordinates) <= 0:
        raise ValueError("COORDINATES must be nonempty")

    for index, granules in enumerate(coordinates):
        if len(granules) <= 0:
            raise ValueError("GRANULES at %d coordinate must be nonempty" % index)

    number_of_variables = len(coordinates)
    probabilities = prepare_probabilities(coordinates[0], number_of_variables)

    discrete_variables = iidfv_discrete(probabilities, number_of_variables)

    variables = [continuous_variable(variable, coordinate)
        for variable, coordinate in izip(discrete_variables, coordinates)]

    return np.array(variables)

