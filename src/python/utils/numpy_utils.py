__all__ = ["as_column", "as_row"]

# creates 2d column (Nx1) from 1d array without copying
def as_column(a):
    if a.ndim != 1:
        raise ValueError("Needed 1d array, got " + "x".join([str(n) for n in a.shape]) + " array")
    return a.reshape((len(a), 1))

# creates 2d row (1xN) from 1d array without copying
def as_row(a):
    if a.ndim != 1:
        raise ValueError("Needed 1d array, got " + "x".join([str(n) for n in a.shape]) + " array")
    return a.reshape((1, len(a)))

