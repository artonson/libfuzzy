from numpy import dot

from numpy_utils import as_column, as_row

def projection_operator(elements):
    """Compute projection to a linear subspace
    spanned over given elements"""
    return dot(as_column(elements), as_row(elements))

