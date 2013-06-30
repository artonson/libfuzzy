
import sys
import numpy as np

__all__ = [
    "prepare_probabilities",
    "random_from_distribution"
]

MAX_FLOAT = int(sys.float_info.max)

# number is long int which can exceed MAX_FLOAT
# multiplier is float
# It is supposed that number * multiplier
# is less than MAX_FLOAT (otherwise an exception
# is raised)
def smart_multiply(base, power, multiplier):
    decompositions = []

    number = base ** power
    while number > MAX_FLOAT:
        if 0 == power % 2:
            decompositions.append(power / 2)
            power /= 2
            number = base ** power
        else:
            decompositions.append(1)
            power -= 1
            number = base ** power

    result = number * multiplier
    while decompositions:
        result *= decompositions.pop()

    return result

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
            multiplier = float(probabilities[current_event - 2])
            numerator -= (smart_multiply(current_event - 1, number_of_variables, multiplier) - \
                smart_multiply(current_event - 2, number_of_variables, multiplier))
        try:
            denominator = float(current_event ** number_of_variables - \
                (current_event - 1) ** number_of_variables)
        except OverflowError:
            probabilities[current_event - 1] = 0.
        else:
            lhs_value = numerator / (denominator + 1.)
            rhs_value = numerator / denominator

            probabilities[current_event - 1] = lhs_value + (rhs_value - lhs_value) / 2.

    if number_of_events > 1:
        multiplier = float(probabilities[current_event - 2])
        numerator -= (smart_multiply(number_of_events - 1, number_of_variables, multiplier) - \
            smart_multiply(number_of_events - 2, number_of_variables, multiplier))
    try:
        denominator = float(number_of_events ** number_of_variables - \
            (number_of_events - 1) ** number_of_variables)
    except OverflowError:
        probabilities[-1] = 0.
    else:
        probabilities[-1] = numerator / denominator

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

