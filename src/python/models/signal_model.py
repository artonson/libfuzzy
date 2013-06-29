from itertools import izip

import numpy as np
from scipy.stats import norm

from fuzzy_random.iidfv import iidfv

__all__ = ["ProbabilisticSignalModel", "FuzzySignalModel",
    "coordinated_fuzzy_model"]

"""
Describe models that are elements
of some n-dimensional space.
"""

class SignalModel(object):
    """
    Arbitrary signal model.
    """
    def __init__(self, *args, **kwargs):
        self.value = kwargs.get("value", None)
        self.coordinated_model = kwargs.get("coordinated_model", None)

    def sample(self):
        """Produces a single realization of a model"""
        raise NotImplementedError()

class ProbabilisticSignalModel(SignalModel):
    """
    An n-dimentional random vector with given
    expectation, covariance matrix (optional),
    and distribution (optional).
    """
    def __init__(self, *args, **kwargs):
        super(ProbabilisticSignalModel, self).__init__(*args, **kwargs)
        self.expectation = kwargs.get("expectation", None)
        self.cov_matrix = kwargs.get("cov_matrix", None)
        self.distribution = kwargs.get("distribution", None)

    def sample(self):
        """Sample a random vector with given parameters.
        Note that currently only iidrv are available
        for uniform distribution.
        """
        if None != self.cov_matrix:
            if None == self.distribution or self.distribution == "normal":
                self.value = np.random.multivariate_normal(self.expectation, self.cov_matrix)
            elif self.distribution == "uniform":
                variances = np.diag(self.cov_matrix)
                self.value = self.expectation + np.random.uniform()
            else:
                raise NotImplementedError("Unknown distribution")
        else:
            self.value = self.expectation

class Granule:
    """
    A granule is a subset of n-dimensional
    space with an assigned probability.
    It is used to designate a simple event
    in a reduced probability space.
    """
    def __init__(self, intervals, probability):
        self.intervals = intervals
        self.probability = probability

    def left_most(self):
        return self.intervals[0][0]

    def right_most(self):
        return self.intervals[-1][1]

class FuzzySignalModel(SignalModel):
    """
    An n-dimentional fuzzy vector with given distribution.
    """
    def __init__(self, *args, **kwargs):
        super(FuzzySignalModel, self).__init__(*args, **kwargs)
        self.expectation = kwargs.get("expectation", None)
        self.coordinates = kwargs.get("coordinates", None)

    def sample(self):
        if None != self.coordinates:
            self.value = iidfv(self.coordinates)
        else:
            self.value = self.expectation


def granulate_normal_pdf(mean, std, ordering, sharpness):
    """Create granules for normal density"""
    if sharpness <= 1.0:
        raise ValueError("sharpness parameter should be greater then 1.0")

    if len(ordering) <= 0:
        raise ValueError("ordering parameter should have nonzero length")

    granules = []
    total_probability = 0.
    for index, relation in enumerate(ordering):
        index += 1
        probability = sharpness ** index / (1. + sharpness) ** index
        total_probability += probability
        cdf_probabilities = [(1. - total_probability) / 2., (1. + total_probability) / 2.]

        intervals = []
        inner = norm.ppf(cdf_probabilities, loc=mean, scale=std)
        if index == 1:
            intervals.append(inner.tolist())
        elif index == len(ordering):
            intervals.append([-1e4, granules[-1].left_most()])
            intervals.append([granules[-1].left_most(), 1e4])
            probability = 1. / (1. + sharpness) ** index
        else:
            intervals.append([inner[0], granules[-1].left_most()])
            intervals.append([granules[-1].right_most(), inner[1]])

        granules.append(Granule(intervals, probability))

    for index in xrange(len(granules)):
        granules[index].probability /= total_probability

    return granules

def coordinated_fuzzy_model(probabilistic_model, **kwargs):
    """Create a fuzzy model maximally coordinated
    with a given probabilistic model.
    coordinated_fuzzy_model granulates the density
    of the supplied probability distribution
    to create a granulated fuzzy model.

    kwargs may contain the following parameters
    - sharpness (number, must be greater 1): determines how
      closely to approximate density (1.0 == maximally)
    - ordering (list, must contain only 1 or 0):
      determines which ordering of possibilities should
      be in the resulting fuzzy distribution. Note that
      as many granules are created, as there are
      items in the ordering list.
    """

    sharpness = kwargs.get("sharpness", 1.1);
    ordering = kwargs.get("ordering", [1, 1, 1])

    if sharpness <= 1.0:
        raise ValueError("sharpness parameter should be greater then 1.0")

    coordinates = None
    if None != probabilistic_model.cov_matrix:
        if None == probabilistic_model.distribution or probabilistic_model.distribution == "normal":
            coordinates = []
            for mean, var in izip(probabilistic_model.expectation, np.diag(probabilistic_model.cov_matrix)):
                coordinates.append(granulate_normal_pdf(mean, np.sqrt(var), ordering, sharpness))
        elif probabilistic_model.distribution == "uniform":
            raise NotImplementedError()
        else:
            raise ValueError("unknown distribution")

    return FuzzySignalModel(probabilistic_model.expectation,
        coordinates=coordinates, coordinated_model=probabilistic_model)

