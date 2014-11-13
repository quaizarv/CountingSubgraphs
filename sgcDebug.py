import struct
import operator
import snap
from collections import *
import random
import numpy as np
from genPrime import *
from kwise import *

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
         
SEL_BITS = 4
class HashFuncSelector:
  # H is the template Graph (aka subgraph) frequency of which is being
  # counted in an input graph G and N is domain of node ids in G
  def __init__(self, H, N):
    self.t = H.GetNodes()
    self.XcHashFnTable = {}
    for c in H.Nodes():
      outDeg = c.GetOutDeg()
      if outDeg not in self.XcHashFnTable:
        thf = TabulatedHashFamily(4*H.GetEdges(), N, outDeg)
        self.XcHashFnTable[outDeg] = [None] * (2 ** SEL_BITS)
        for i in range((2 ** SEL_BITS)):
          self.XcHashFnTable[outDeg][i] = thf.getHashFunc()

    self.YHashFnTable = [None] * (2 ** SEL_BITS)
    thf = TabulatedHashFamily(4*H.GetEdges(), N, H.GetNodes())
    for i in range((2 ** SEL_BITS)):
      self.YHashFnTable[i] = thf.getHashFunc()

  def getXc(self, c, selState):
      cId = c.GetId()
      outDeg = c.GetOutDeg()
      print "getXc: ", outDeg, selState[cId]
      return self.XcHashFnTable[outDeg][selState[cId]]

  def getY(self, selState):
      print "getY: ", selState[self.t]
      return self.YHashFnTable[selState[self.t]]

  def getQ(self, selState):
    return selState[self.t+1]

  def selectionStateCreate(self):
    bits = random.randint(0, 2**32 - 1) / 2 ** (32 - (SEL_BITS*(self.t+1) +
                                                      self.t))
    selector = {}
    # last t bits are used for selecting Q
    selector[self.t+1] = bits % (2**self.t)
    print "bits, Q:", bits, selector[self.t+1]
    bits = bits / (2**self.t)

    # (4t+1 to 4t+4) bits are used for selecting hash function for Y
    selector[self.t] = bits % (2 ** SEL_BITS)
    bits = bits / (2 ** SEL_BITS)

    # first t bits are used in selecting hash functions for Y
    for c in range(self.t-1,-1, -1):
      selector[c] = bits % (2 ** SEL_BITS)
      bits = bits / (2 ** SEL_BITS)
    return selector
      
    
class Estimator:
  # H is the template Graph (aka subgraph) frequency of which is being
  # counted in the input graph G
  def __init__(self, H, hashFnSelector):
    self.hashFnSelState = hashFnSelector.selectionStateCreate()

    def createXc(c, XcHash):
      ## precompute a few things so that we don't have to compute them often
      outDeg = c.GetOutDeg()
      roots = unityRoots(outDeg)
      def Xc(u):
        return roots[XcHash(u)]
      return Xc  ## return the closure

    def computeXcTable(H):
      XcTable = {}
      for c in H.Nodes():
        cId = c.GetId()
        XcHash = hashFnSelector.getXc(c, self.hashFnSelState)
        XcTable[cId] = createXc(c, XcHash)
      return XcTable    

    def createY(H):
      ## precompute a few things so that we don't have to compute them often
      YHashFn = hashFnSelector.getY(self.hashFnSelState)
      def Y(u):
        return 2 ** (YHashFn(u))
      return Y

    self.H = H
    tau = 2**H.GetNodes() - 1
    self.Q = unityRoots(tau)[hashFnSelector.getQ(self.hashFnSelState) % tau]
    print self.Q
    self.XcTable = computeXcTable(H)
    self.Y = createY(H)
    self.Z = defaultdict(lambda: 0)
    self.degHNode = {n.GetId():n.GetOutDeg() for n in H.Nodes()}

  # 'edge' is the input edge from G
  def updateZ(self, gEdge):
    def matchProb((a, b), (u, v)):
      Xa = self.XcTable[a]
      Xb = self.XcTable[b]
      return (Xa(u) * Xb(v) * 
              (self.Q ** (1.0 * self.Y(u)/self.degHNode[a])) *
              (self.Q ** (1.0 * self.Y(v)/self.degHNode[b])))
      
    u, v = gEdge.GetId()
    for hEdge in self.H.Edges():
      a, b = hEdge.GetId()
      self.Z[(a, b)] += matchProb((a,b), (u,v)) + matchProb((a, b), (v, u))

  def clearZ(self):
    for hEdge in self.H.Edges():
      a, b = hEdge.GetId()
      self.Z[(a, b)] = 0

  def debug(self, G):
    def matchProb((a, b), (u, v)):
      Xa = self.XcTable[a]
      Xb = self.XcTable[b]
      #print "Xa Xb: ", (a, b), (u, v), Xa(u), Xb(v)
      return (Xa(u) * Xb(v) * 
              (self.Q ** (1.0 * self.Y(u)/self.degHNode[a])) *
              (self.Q ** (1.0 * self.Y(v)/self.degHNode[b])))

    def edgeOrientations(e):
      u, v = e
      return [(u,v), (v, u)]

    def edgeSetOrientations(T):
      return [(o1, o2, o3)
              for o1 in edgeOrientations(T[0])
              for o2 in edgeOrientations(T[1])
              for o3 in edgeOrientations(T[2])]

    Z = 0
    edges = [e.GetId() for e in G.Edges()]
    hEdges = tuple([e.GetId() for e in self.H.Edges()])
    m = len(edges)
    for i in range(m):
      for j in range(m):
        for k in range(m):
          T = [edges[ix] for ix in [i,j,k]]
          #print T
          for tEdgeSet in edgeSetOrientations(T):
            prod = 1
            for ((a,b), (u,v)) in zip(hEdges, tEdgeSet):
              M = matchProb((a, b), (u, v))
              #print (a,b), (u, v), M
              prod = prod * M
            #print "Z for edgeSet: ", tEdgeSet,
            #print "is: ", prod
            Z += prod
    #print "Total: ", Z
    return Z

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

def countTriangles(G):
  H = snap.GenFull(snap.PUNGraph, 3)
  hashFnSel = HashFuncSelector(H, 2**16)
  total = 0.0
  results = []
  estCount = 1
  for i in range(estCount):
    est = Estimator(H, hashFnSel)
    for e in G.Edges():
      est.updateZ(e)
    results.append(est.subgraphCount())
    print i, results[-1]
  return results, sum(results)/estCount


def debugTriangles(G):
  H = snap.GenFull(snap.PUNGraph, 3)
  hashFnSel = HashFuncSelector(H, 2**4)
  total = 0.0
  results = []
  estCount = 100
  for i in range(estCount):
    est = Estimator(H, hashFnSel)
    for e in G.Edges():
      est.updateZ(e)
    #results.append(est.debug(G))
    results.append(est.subgraphCount())
    print i, results[-1]
  return results, sum(results)/estCount


