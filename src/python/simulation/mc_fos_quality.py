#!/usr/bin/python
import numpy as np

from itertools import product
import os
import sys

sys.path[1:1] = [os.path.normpath(os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    ".."))]

#@mr_include_dirs
from measurement.scheme import LinearMeasurementScheme
from models.signal_model import ProbabilisticSignalModel, FuzzySignalModel, coordinated_fuzzy_model
from reduction.reduction import ProbabilisticLinearReduction, FuzzyLinearReduction
from signals.example_signals import fos_eigenvector
from sensors.first_order_sensor import first_order_sensor
from utils.linalg_utils import projection_operator
#@end

def probabilistic_quality():
    DIMENSIONS = 4
    signal_expectation = fos_eigenvector(1, DIMENSIONS)
    signal_model = ProbabilisticSignalModel(expectation=signal_expectation)

    NOISE_STD = 10.0
    noise_cov_matrix = NOISE_STD * np.eye(DIMENSIONS)
    noise_model = ProbabilisticSignalModel(
        expectation=np.zeros(signal_expectation.shape),
        cov_matrix=noise_cov_matrix)

    alpha = 0.1
    beta = 0.2
    operator = first_order_sensor(alpha, beta, DIMENSIONS)
    probabilistic_scheme = LinearMeasurementScheme(operator, signal_model, noise_model)
    probabilistic_scheme.measure()

    ideal_operator = projection_operator(signal_expectation) / DIMENSIONS
    reduction = ProbabilisticLinearReduction(probabilistic_scheme, ideal_operator)
    reduction.compute()
    print reduction.error.value

def fuzzy_quality():
    DIMENSIONS = 4
    signal_expectation = fos_eigenvector(1, DIMENSIONS)
    signal_model = FuzzySignalModel(expectation=signal_expectation)

    NOISE_STD = 10.0
    noise_cov_matrix = NOISE_STD * np.eye(DIMENSIONS)
    noise_model = ProbabilisticSignalModel(
        expectation=np.zeros(signal_expectation.shape),
        cov_matrix=noise_cov_matrix)

    fuzzy_noise_model = coordinated_fuzzy_model(noise_model)

    alphas = np.linspace(-1., +1., 100)
    betas = np.linspace(-1., +1., 100)
    for alpha, beta in product(alphas, betas):
        operator = first_order_sensor(alpha, beta, DIMENSIONS)
        fuzzy_scheme = LinearMeasurementScheme(operator, signal_model, fuzzy_noise_model)
        fuzzy_scheme.measure()

        ideal_operator = projection_operator(signal_expectation) / DIMENSIONS
        reduction = FuzzyLinearReduction(fuzzy_scheme, ideal_operator)
        reduction.compute()
        print reduction.error.value

def main():
    probabilistic_quality()
    fuzzy_quality()

if __name__ == "__main__":
    main()

