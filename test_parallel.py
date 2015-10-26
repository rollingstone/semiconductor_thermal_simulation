__author__ = 'kamal'



import sys
import os

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from multiprocessing.pool import Pool
# import distribute
from multiprocessing.process import Process
from multiprocessing.queues import Queue

# from multiprocessing.Process import pro



v1 = np.random.rand(4, 4)
v2 = np.random.rand(4, 4)



def func(v1):
    return v1 * v1


def mmul(matrix):
    for i in range(10):
        matrix = matrix * matrix
    return matrix


matrices = []

for i in range(4):
    matrices.append(np.random.random_integers(100, size=(1000, 1000)))



p = Pool(processes=4)

p.apply_async(mmul, matrices)
p.join()
p.close()



