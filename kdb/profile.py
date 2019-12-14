import array
import math
import functools

import logging
logger = logging.getLogger(__file__)

def magnitude(x):
    """ Calculates the euclidean norm of the unit vector of the input

    :param x: The k-mer profile array (spectrum)
    :type x: array.array
    :returns: 
    :rtype: float

    """
    if not isinstance(x, array.array):
        raise TypeError("kdb.profile.magnitude expects an array as its positional argument")
    return math.sqrt(functools.reduce(lambda a,b: a+(b**2), x, 0))

def euclidean_distance(x, y):
    """The euclidean distance of the unit vectors of x and y

    :param x: 
    :type x: array.array
    :param y: 
    :type y: array.array
    :returns: 
    :rtype: float

    """
    if not isinstance(x, array.array):
        raise TypeError("kdb.profile.euclidean_distance expects an array as its first positional argument")
    elif not isinstance(y, array.array):
        raise TypeError("kdb.profile.euclidean_distance expects an array as its second positional argument")
    elif len(x) != len(y):
        raise TypeError("kdb.profile.euclidean_distance expects arrays of equal lengths")
    magnitude_x = magnitude(x)
    magnitude_y = magnitude(y)
    z = 0
    for i in range(len(x)):
        z += (x[i]/magnitude_x - y[i]/magnitude_y)**2
    return math.sqrt(z)

def correlation_distance(x, y):
    """The correlation distance of the unit vectors of x and y

    :param x: 
    :type x: array.array
    :param y: 
    :type y: array.array
    :returns: 
    :rtype: float

    """
    if not isinstance(x, array.array):
        raise TypeError("kdb.profile.correlation_distance expects an array as its first positional argument")
    elif not isinstance(y, array.array):
        raise TypeError("kdb.profile.correlation_distance expects an array as its second positional argument")
    elif len(x) != len(y):
        raise TypeError("kdb.profile.correlation_distance expects arrays of equal lengths")
    n = len(x)
    magnitude_x = magnitude(x)
    magnitude_y = magnitude(y)
    mean_x = functools.reduce(lambda a,b: a+(b/magnitude_x), x, 0) / n
    mean_y = functools.reduce(lambda a,b: a+(b/magnitude_y), y, 0) / n

    ssxx = functools.reduce(lambda a,b: a+((b/magnitude_x) - mean_x)**2, x, 0)
    # ssxx = math.sqrt(functools.reduce(lambda a,b: a+(b/magnitude_x)**2, x)) - (n * mean_x * mean_x)
    ssyy = functools.reduce(lambda a,b: a+((b/magnitude_y) - mean_y)**2, y, 0)
    # ssyy = math.sqrt(functools.reduce(lambda a,b: a+(b/magnitude_y)**2, y)) - (n * mean_y * mean_y)
    ssxy = 0
    for i in range(n):
        ssxy += ((x[i]/magnitude_x - mean_x) * (y[i]/magnitude_y - mean_y))
    return ssxy/math.sqrt(ssxx*ssyy)
    

