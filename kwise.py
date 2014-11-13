from collections import *
import numpy as np
import random
import matplotlib.pyplot as plt

def hashValue(poly, x, p, m):
  return np.polyval(poly, x) % p % m

def hashRangeProb(k, hashpoly = None, p=None, domainSize=None, rangeSize=None):
  if p is None: p = 10007
  if domainSize is None: domainSize = p
  if rangeSize is None: rangeSize = p
  if hashpoly is None:
    hashpoly = randomNDegreePoly(k-1, p)
  return Counter([hashValue(hashpoly, x, p, rangeSize) 
                  for x in range(domainSize)])

def plotCounter(cntr):
  plt.plot(cntr.keys(), cntr.values())
  plt.show()


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
  p = 5
  #p = 10007  ## Prime >> m
  d = defaultdict(lambda: 0)
  hashpoly = randomNDegreePoly(k-1, p)
  for t in generateKTuples(k, range(domainSize)):
    if not tupleElementsUnique(t): continue
    d[tuple([hashValue(hashpoly, x, p, rangeSize) for x in t])] += 1
    if sampleSize and sum(d.values()) == sampleSize: break
  return d

def test():
  testKWiseIndependence(3, 64, 2, None)
  #testKWiseIndependence(4, 32, 2, None)

########################################################################

def mChooseN(m, n):
  if m <= n or m is 0 or n is 0: return 1
  else: return mChooseN(m-1, n-1) + mChooseN(m-1, n)

# Kwise Hash based on the paper by Martin Dietzfelbinger - Universal Hashing and
# k-wise independent random variables via integer arithmetic without primes
class MDHashFamily:
  def __init__(self, k, domainSize, rangeSize):
    self.k = k
    self.u = domainSize
    self.m = rangeSize
    self.l = (self.u-1)**mChooseN(self.k,2)
    self.r = self.m * self.l

  def getHashFunc(self):
    coeffs = [random.randint(0, self.r) for i in  range(self.k-1)]
    coeffs.append(random.randint(1, self.r))
    return lambda x: (sum([(coeffs[i] * (x ** i)) % self.r
                           for i in range(self.k)]) % self.r) / self.l

class TabulatedHashFamily:
  def __init__(self, k, domainSize, rangeSize):
    if domainSize <= 2**16: q = 1
    else: q = 2
    self.k = k
    self.q = q
    self.u = domainSize
    self.m = rangeSize
    self.hashFnsCount = (k-1)*(q-1) + 1
    r = (k-2)*(q-1)

    if self.q == 1:
      self.H = MDHashFamily(self.k, self.u, self.m)
    else:
      self.H = MDHashFamily(self.k, 2 ** 16, self.m)
    self.V = np.transpose(np.vander(range(q+r), q, True))

  def getHashFunc(self):
    tables = [None] * self.hashFnsCount
    for i in range(self.hashFnsCount):
      h = self.H.getHashFunc()
      if self.q == 1:
        tables[i] = [h(x) for x in range(self.u)]
      else:
        tables[i] = [h(x) for x in range(2**16)]

    def hashFunc(x):
      # Shortcut if the key is only one character long
      if self.hashFnsCount == 1: return tables[0][x]

      # Multi-character key
      x_v = [x / (2**16), x % (2**16)]
      z_v = np.dot(x_v, self.V) % (2**16)
      val = 0
      for i in range(self.hashFnsCount):
        val = val ^ tables[i][z_v[i]]
      return val
      
    return hashFunc
