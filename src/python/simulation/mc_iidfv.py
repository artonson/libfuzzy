#!/usr/bin/python
import numpy as np

from itertools import product, repeat, izip
import multiprocessing as mp
import os
import sys
import time

sys.path[1:1] = [os.path.normpath(os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    ".."))]

from models.signal_model import ProbabilisticSignalModel, FuzzySignalModel, coordinated_fuzzy_model

VARS_COUNT = np.arange(2,1000)
MC_ITERATIONS = int(1e3)

def time_it(f, *args):
    start = time.clock()
    f(*args)
    return time.clock() - start

def probabilistic_iidfv(dimensions):
    import numpy as np

    NOISE_STD = 1.0
    signal_expectation = np.zeros((dimensions,))
    noise_cov_matrix = NOISE_STD * np.eye(dimensions)
    noise_model = ProbabilisticSignalModel(
        expectation=np.zeros(signal_expectation.shape),
        cov_matrix=noise_cov_matrix)

    times = []
    for iteration in xrange(MC_ITERATIONS):
        times.append(time_it(lambda model: model.sample(), noise_model))

    return dimensions, np.mean(times), np.std(times)

def fuzzy_iidfv(dimensions):
    import numpy as np

    NOISE_STD = 1.0
    signal_expectation = np.zeros((dimensions,))
    noise_cov_matrix = NOISE_STD * np.eye(dimensions)
    noise_model = ProbabilisticSignalModel(
        expectation=np.zeros(signal_expectation.shape),
        cov_matrix=noise_cov_matrix)
    fuzzy_noise_model = coordinated_fuzzy_model(noise_model)

    times = []
    for iteration in xrange(MC_ITERATIONS):
        times.append(time_it(lambda model: model.sample(), fuzzy_noise_model))

    return dimensions, np.mean(times), np.std(times)

def do_test(func):
    pool = mp.Pool(processes=24)
    for result in pool.imap_unordered(func, VARS_COUNT):
        print "\t".join([str(value) for value in result])

if __name__ == "__main__":
    if len(sys.argv) < 2 or len(sys.argv) == 2 and sys.argv[1] == '-h':
        print "Usage: ./mc_iidfv.py [prob | poss]"
        sys.exit(1)

    test = sys.argv[1]
    if test == "prob":
        do_test(probabilistic_iidfv)
    elif test == "poss":
        do_test(fuzzy_iidfv)
    else:
        raise ValueError("Unknown test: %s", test)

