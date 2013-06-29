
import numpy as np

__all__ = [
    "prepare_probabilities",
    "random_from_distribution"
]

def prepare_probabilities(ordering, number_of_variables):
    if number_of_variables <= 0:
        raise ValueError("NUMBER_OF_VARIABLES must be a positive number")

    if len(ordering) <= 0:
        raise ValueError("PROBABILITIES must be nonempty")

    number_of_events = len(ordering)
    probabilities = np.zeros(shape=(number_of_events,))

    numerator = 1.
    for current_event in xrange(1, number_of_events):
        if current_event > 1:
            numerator -= float(((current_event - 1) ** number_of_variables - \
                (current_event - 2) ** number_of_variables)) * \
                probabilities[current_event - 2]
        denominator = current_event ** number_of_variables - \
            (current_event - 1) ** number_of_variables

        lhs_value = float(numerator) / float(denominator + 1.)
        rhs_value = float(numerator) / float(denominator)

        probabilities[current_event - 1] = lhs_value + (rhs_value - lhs_value) / 2.

    if number_of_events > 1:
        numerator -= float(((number_of_events - 1) ** number_of_variables - \
            (number_of_events - 2) ** number_of_variables)) * \
            probabilities[number_of_events - 2]
    denominator = number_of_events ** number_of_variables - \
        (number_of_events - 1) ** number_of_variables
    probabilities[-1] = float(numerator) / float(denominator)

    return probabilities

def random_from_distribution(distribution):
    value = np.random.random()
#    print"np.random.random:", value
    sum_value = 0
    for index, probability in enumerate(distribution):
        sum_value += probability
        if value <= sum_value:
            return index
    return index

