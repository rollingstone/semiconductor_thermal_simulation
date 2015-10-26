__author__ = 'kamal'


import numpy as np
import scipy as sci
import scipy.linalg as la
import cProfile, pstats, StringIO
import re



def func1(val):
    av = np.random.randint(0,100,(50,50))
    vv = []
    for i in xrange(val):
        vv = la.inv(av)

    return vv




pr = cProfile.Profile()
pr.enable()
vv = func1(10000)
pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print s.getvalue()


# cProfile.run('re.compile("func1")')