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
    selector = {}
    # Q value
    selector[self.t+1] = random.randint(0, 2**self.t - 1)
    print "Q:", selector[self.t+1]

    # Y
    selector[self.t] = random.randint(0, 2**SEL_BITS - 1)

    # Xc
    for c in range(self.t):
      selector[c] = random.randint(0, 2**SEL_BITS - 1)
    return selector
      
    
class Estimator:
  # H is the template Graph (aka subgraph) frequency of which is being
  # counted in the input graph G
  def __init__(self, G, H, hashFnSelector):
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

    def createMTable():
      assert(self.G and self.H and self.XcTable and self.Q and self.Y and
             self.degHNode)
      MTable = {}
      for u in range(self.G.GetNodes()):
        for c in range(self.H.GetNodes()):
          Xc = self.XcTable[c]
          MTable[(u,c)] = Xc(u) * (self.Q ** (1.0 * self.Y(u)/self.degHNode[c]))
      return MTable

    def powerOfQ(u, c):
      return (self.Q ** (1.0 * self.Y(u)/self.degHNode[c]))
      
    self.G = G
    self.H = H
    tau = 2**H.GetNodes() - 1
    self.Q = unityRoots(tau)[hashFnSelector.getQ(self.hashFnSelState) % tau]
    print self.Q
    self.XcTable = computeXcTable(H)
    self.Y = createY(H)
    self.Z = defaultdict(lambda: 0)
    self.degHNode = {n.GetId():n.GetOutDeg() for n in H.Nodes()}
    self.MTable = createMTable()

  # 'edge' is the input edge from G
  def updateZ(self, gEdge):
    def matchProb((a, b), (u, v)):
      return self.MTable[(u, a)] * self.MTable[(v, b)]
      Xa = self.XcTable[a]
      Xb = self.XcTable[b]
      return Xa(u) * Xb(v) * \
          (self.Q ** (1.0 * self.Y(u)/self.degHNode[a])) * \
          (self.Q ** (1.0 * self.Y(v)/self.degHNode[b]))

    u, v = gEdge.GetId()
    for hEdge in self.H.Edges():
      a, b = hEdge.GetId()
      self.Z[(a, b)] += matchProb((a,b), (u,v)) + matchProb((a, b), (v, u))

  def clearZ(self):
    for hEdge in self.H.Edges():
      a, b = hEdge.GetId()
      self.Z[(a, b)] = 0

  """    
  def edgeOrientations(e):
    u, v = e
    return [(u,v), (v, u)]

  def edgeSetOrientations(T):
    return [(o1, o2, o3)
            for o1 in edgeOrientations(T[0])
            for o2 in edgeOrientations(T[0])
            for o3 in edgeOrientations(T[0])]

  def debug(self, G):
    edges = [e.GetId() for e in G.Edges()]
    m = len(edges)
    for i in range(m):
      for j in range(m):
        for k in range(m):
          T = [edges[ix] for ix in [i,j,k]]
          for o1, o2, o3 in edgeSetOrientations(T):
            
  """
          
          

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
  N = 2**16
  assert(G.GetNodes() <= N)
  hashFnSel = HashFuncSelector(H, N)
  total = 0.0
  results = []
  estCount = 3000
  for i in range(estCount):
    est = Estimator(G, H, hashFnSel)
    for e in G.Edges():
      est.updateZ(e)
    results.append(est.subgraphCount())
    print i, results[-1]
  return int(sum(results)/estCount), results


