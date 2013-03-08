#!/usr/bin/env python
#from math import sqrt, ceil
import numpy as np

def rwh_pcn(n):
    # http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    """ Input n>=6, Returns a list of primes, 2 <= p < n """
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    for i in xrange(1,int(n**0.5)/3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[       k*k/3     ::2*k] = False
            sieve[k*(k-2*(i&1)+4)/3::2*k] = False
    return 1 + np.count_nonzero(sieve)
    #return np.r_[2,3,((3*np.nonzero(sieve)[0]+1)|1)]


print rwh_pcn(800000000)
