import math


def array_almost_equals(a, b, threshold=1e-3):
    """Check if two arrays containing float number are almost equals"""

    if len(a) != len(b):
        raise RuntimeError()

    return all(math.fabs(a[i] - b[i]) < threshold for i in range(len(a)))
