import struct
import operator
import snap
from collections import *
import random
import numpy as np
import mmh3
from genPrime import *

def permutations(l):
  if l == []: return [[]]
  else:
    return [l[i:i+1] + p 
            for i in range(len(l))
            for p in permutations(l[0:i] + l[i+1:len(l)])]

def test_permutations():
  l = [1, 2, 3]
  assert(tuple(permutations(l)) ==
         tuple([[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3, 1, 2],
                [3, 2, 1]]))

def isAuto(H, perm):
  N = H.GetNodes()
  V = [n.GetId() for n in H.Nodes()]
  for i in range(N):
    for j in range(i+1, N):
      if H.IsEdge(V[i], V[j]) and not H.IsEdge(perm[i], perm[j]):
        return False
  return True

def automorphismCount(H):
  N = H.GetNodes()
  return sum([1 for p in permutations(range(N)) if isAuto(H, p)])

def test_automorphismCount():
  H = snap.GenGrid(snap.PUNGraph, 3, 2)
  assert(automorphismCount(H) == 4)
  H = snap.GenStar(snap.PUNGraph, 5)
  assert(automorphismCount(H) == 24)

def unityRoots(n):
  coeff = [0 for i in range(n+1)]
  coeff[0] = 1
  coeff[n] = -1
  return np.roots(coeff)

EPSILON = 0.0000001
def test_unityRoots():
  assert((unityRoots(2) == np.array([-1.,  1.])).all())
  assert(abs(unityRoots(3)[0].real + 0.5) < EPSILON)
  assert(abs(unityRoots(3)[0].imag - sqrt(3)/2) < EPSILON)
         
def kWiseHash(n, seed):
  return mmh3.hash(struct.pack('I', n), seed[0])

def kWiseHashSeed(k, p):
  return [random.randint(1, p) for _ in range(k)]

def createXc(c, H, p):
  ## precompute a few things so that we don't have to compute them often
  seed = kWiseHashSeed(4 * H.GetEdges(), p)
  roots = unityRoots(c.GetOutDeg())
  outDeg = c.GetOutDeg()
  def Xc(u):
    return roots[kWiseHash(u, seed) % outDeg]
  return Xc  ## return the closure

def computeXcTable(H, N):
  # Here nId is a node Id in template graph and N is the # of nodes
  # in the input Graph G
  p = firstPrimeAfter(N)
  XcTable = {}
  for c in H.Nodes():
    cId = c.GetId()
    XcTable[cId] = createXc(c, H, p)
  return XcTable    

def createY(H, N):
  ## precompute a few things so that we don't have to compute them often
  p = firstPrimeAfter(N)
  seed = kWiseHashSeed(4 * H.GetEdges(), p)
  t = H.GetNodes()
  def Y(u):
    return 2 ** (kWiseHash(u, seed) % t)
  return Y ## return the closure

class Estimator:
  # H is the template Graph (aka subgraph) frequency of which is being
  # counted in the input graph G and N is size of G
  def __init__(self, H, N):
    self.H = H
    tau = 2**H.GetNodes() - 1
    self.Q = unityRoots(tau)[random.randint(0, tau-1)]
    self.XcTable = computeXcTable(H, N)
    self.Y = createY(H, N)
    self.Z = defaultdict(lambda: 0)
    self.degHNode = {n.GetId():n.GetOutDeg() for n in H.Nodes()}

  # 'edge' is the input edge from G
  def updateZ(self, gEdge):
    def matchProb((a, b), (u, v)):
      Xa = self.XcTable[a]
      Xb = self.XcTable[b]
      return (Xa(u) * Xb(v) * 
              (self.Q ** (self.Y(u)/self.degHNode[a])) *
              (self.Q ** (self.Y(v)/self.degHNode[b])))
      
    u, v = gEdge.GetId()
    for hEdge in self.H.Edges():
      a, b = hEdge.GetId()
      self.Z[(a, b)] += matchProb((a,b), (u,v)) + matchProb((a, b), (v, u))

  def clearZ(self):
    for hEdge in self.H.Edges():
      a, b = hEdge.GetId()
      self.Z[(a, b)] = 0
    
  def debug(self, G):
    edges = [e.GetId() for e in G.Edges()]
    m = len(edges)
    for i in range(m):
      for j in range(m):
        for k in range(m):
          e1 = edges[i]
          e2 = edges[j]
          e3 = edges[k]
          
          

  def subgraphCount(self):
    t = self.H.GetNodes()
    result = (reduce(operator.mul, [t * 1.0/i for i in range(1, t+1)])/
              automorphismCount(self.H))
    result = reduce(operator.mul,
                    [self.Z[e.GetId()] for e in self.H.Edges()], result)
    return result.real

def countTrianglesBruteForce(G):
  count = 0
  for i in range(G.GetNodes()):
    for j in range(i+1, G.GetNodes()):
      for k in range(j+1, G.GetNodes()):
        if G.IsEdge(i,j) and G.IsEdge(j, k) and G.IsEdge(k, i):
          count +=1
  return count

def countSubgraphs(G):
  total = 0.0
  results = []
  estCount = 1
  for i in range(estCount):
    H = snap.GenRndGnm(snap.PUNGraph, 3, 3)
    est = Estimator(H, G.GetNodes())
    i = 0
    for e in G.Edges():
      i = i + 1
      if i % 10000 == 0: print i
      est.updateZ(e)
    results.append(est.subgraphCount())
  return sum(results)/estCount


