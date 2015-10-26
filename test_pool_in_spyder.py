# -*- coding: utf-8 -*-
"""
Created on Tue Oct 13 17:38:14 2015

@author: kamal
"""

__author__ = 'kamal'



import sys
import os

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from multiprocessing.pool import Pool
from multiprocessing.process import Process
from multiprocessing.queues import Queue

from parutils import distribute
import sharedmem

import concurrent.futures

from timeit import timeit


v1 = np.random.rand(4, 4)
v2 = np.random.rand(4, 4)



def func(v1):
    return v1 * v1


def mmul(matrix):
    for i in range(20):
        matrix = matrix * matrix
    return matrix

matrices = []

for i in range(4):
    matrices.append(np.random.random_integers(100, size=(1000, 1000)))

print timeit(lambda: map(mmul, matrices), number=20)

p = Pool(processes=6)

print timeit(lambda: p.map(mmul, matrices), number=20)
print timeit(lambda: p.map(mmul, matrices), number=20)
print timeit(lambda: p.map(mmul, matrices), number=20)
print timeit(lambda: p.map(mmul, matrices), number=20)
p.close()
p.join()



