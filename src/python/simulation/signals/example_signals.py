
import numpy as np

from utils.numpy_utils import as_column

def fos_eigenvector(index, dimensions):
    grid = np.linspace(1./dimensions, 1.0, dimensions)
    return np.sqrt(2) * np.sin(np.pi * (1 - grid) / 2)

