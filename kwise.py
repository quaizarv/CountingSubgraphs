from collections import *
import numpy as np
import random


def hashValue(poly, x, p, m):
  return np.polyval(poly, x) % p % m

def generateKTuples(k, domain):
  def empty(stack): return len(stack) == 0
  stack = [domain[0] for _ in range(k)]
  while True:
    yield tuple(stack)
    if stack[-1] < domain[-1]:
      stack[-1] += 1
      continue

    while not empty(stack) and stack[-1] >= domain[-1]: stack.pop()
    if empty(stack): break
    stack[-1] += 1
    while len(stack) < k: stack.append(domain[0])

def test_generateKTuples():
  kTups = [t for t in generateKTuples(3, [0, 1])]
  assert(kTups == [(0, 0, 0),(0, 0, 1),(0, 1, 0),(0, 1, 1),
                   (1, 0, 0),(1, 0, 1),(1, 1, 0),(1, 1, 1)])
  kTups = [t for t in generateKTuples(4, [0, 1])]
  assert(len(kTups) == 16)
  assert(kTups[3] == (0, 0, 1, 1))

def tupleElementsUnique(t):
  return len(set(t)) == len(t)

def randomNDegreePoly(n, p):
  """Returns a degree n polynomial with random coefficients chosen
  from range [1, p-1]
  """
  return [random.randint(1, p-1) for _ in range(n+1)]

def testKWiseIndependence(k, domainSize, rangeSize, sampleSize=None):
  p = 10007  ## Prime >> m
  d = defaultdict(lambda: 0)
  hashpoly = randomNDegreePoly(k-1, p)
  for t in generateKTuples(k, range(domainSize)):
    if not tupleElementsUnique(t): continue
    d[tuple([hashValue(hashpoly, x, p, rangeSize) for x in t])] += 1
    if sampleSize and sum(d.values()) == sampleSize: break
  return d

def test():
  #testKWiseIndependence(3, 100, 2, None)
  testKWiseIndependence(4, 32, 2, None)

