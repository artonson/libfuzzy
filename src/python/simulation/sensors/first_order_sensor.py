
import numpy as np

from utils.numpy_utils import as_column, as_row

def first_order_sensor(alpha, beta, dimensions):
    """Creates matrix for the first order sensor given
    its parameters and dimensions
    """
    integration_step = 1. / dimensions;
    grid = as_column(np.linspace(0., 1., dimensions + 1))
    x1 = grid[1 : dimensions + 1]
    x2 = grid[0 : dimensions]

    temp = as_row(np.ones(dimensions))
    t1 = (np.dot(x1, temp) - np.dot(temp.T, x1.T)) * beta / alpha
    t2 = (np.dot(x2, temp) - np.dot(temp.T, x1.T)) * beta / alpha
    t3 = (np.dot(x1, temp) - np.dot(temp.T, x2.T)) * beta / alpha
    t4 = (np.dot(x2, temp) - np.dot(temp.T, x2.T)) * beta / alpha

    CONST = -alpha / beta**2 / integration_step

    d = CONST * (np.exp(-t1) - np.exp(-t2) - np.exp(-t3) + np.exp(-t4))
    di = 1. / beta + CONST * (1. - np.exp(-beta / alpha * integration_step))
    fos = d - np.triu(d) + di * np.eye(dimensions)

    return fos

