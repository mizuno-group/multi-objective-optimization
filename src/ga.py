import random
from deap import base, creator, tools, algorithms
import itertools

def not_dup(a, b, list):
    k = 0
    while k < 1:
        s = random.randint(a,b)
        if s not in list:
            k = 100
    return s

def rand_nodup(a, b, k):
    if abs(a) + b < k:
        raise ValueError
    r = set()
    while len(r) < k:
        r.add(random.randint(a, b))
    return r
