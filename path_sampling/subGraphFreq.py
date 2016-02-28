import numpy as np
import math
from scipy.misc import comb
import snap
import bisect
import random
import os

def countTriads(G):
  tt = 0
  tp = 0
  for n in G.Nodes():
    t = snap.GetNodeTriads(G, n.GetId())
    tt += t
    #c = snap.GetNodeClustCf(G, n.GetId())
    d = n.GetOutDeg()
    tp = tp + (d * (d-1))/2.0 - t
    #nt = nt + ((d * (d-1))/2.0 - t)
    #nt = nt + (1 - c)* (d * (d-1))/2.0
  return int(tt/3), int(tp)


def cumulativeEdgeDist(G):
  w = np.zeros(G.GetEdges())
  edgeMap = {}
  i = 0
  for e in G.Edges():
    u, v = e.GetId()
    if u == v:
      continue
    t = (G.GetNI(u).GetOutDeg() - 1) * (G.GetNI(v).GetOutDeg() - 1)
    w[i] = t
    edgeMap[i] = (u,v)
    i = i + 1
  W = sum(w)
  w1 = np.cumsum(w/W)
  return w1, W, edgeMap


def sampleNbr(nodeNbrs, exceptNode):
  if len(nodeNbrs) == 1:
    return None
  while True:
    i = random.randint(0, len(nodeNbrs)-1)
    u1 = nodeNbrs[i]
    if u1 != exceptNode:
      break
  return u1
  
def sample(G, cumDist, edgeMap, nbrs):
  u, v = edgeMap[bisect.bisect(cumDist, random.random())]
  u1 = sampleNbr(nbrs[u], v)
  v1 = sampleNbr(nbrs[v], u)
  if u1 == v1: return None
  return (u1, u, v, v1)

def inducedMotif(G, S_l):
  u1, u, v, v1 = S_l
  cycle = False
  edgeCount = 3
  if G.IsEdge(u1, v1):
    cycle = True
    edgeCount += 1
  if G.IsEdge(u1, v):
    edgeCount += 1
  if G.IsEdge(u, v1):
    edgeCount += 1
  if edgeCount == 6: return 6
  elif edgeCount == 5: return 5
  elif edgeCount == 3: return 2
  else: # edgecount is 4
    if cycle: return 4
  return 3

A = np.array([[1, 0, 1, 0, 2, 4],
              [0, 1, 2, 4, 6, 12],
              [0, 0, 1, 0, 4, 12],
              [0, 0, 0, 1, 1, 3],
              [0, 0, 0, 0, 1, 6],
              [0, 0, 0, 0, 0, 1]])

def computeNbrs(G):
  nbrs = {}
  for n in G.Nodes():
    nbrs[n.GetId()] = [nbr for nbr in n.GetOutEdges()]
  return nbrs

def threePathSampler(G, k, cumDist, W, edgeMap, nbrs):
  #cumDist, W, edgeMap = cumulativeEdgeDist(G)
  mCount = [0]*6
  cBar = [0]*6
  i = 0
  for i in range(k):
    S_l = sample(G, cumDist, edgeMap, nbrs)
    if S_l == None: continue
    mCount[inducedMotif(G, S_l)-1] += 1
    if i % 1000 == 0: print inducedMotif(G, S_l), i
  for i in range(2, 6+1):
    cBar[i-1] = (1.0 * mCount[i-1]/k) * W/A[2-1, i-1]
  N_1 = sum([comb(n.GetOutDeg(), 3) for n in G.Nodes()])
  cBar[1-1] = N_1 - cBar[3-1] - 2*cBar[5-1] - 4*cBar[6-1]
  return cBar

def sampleHigherDegNbr(sortedAdjList, v):
  nbrMap, nbrs = sortedAdjList
  idx = nbrMap[v]
  i = random.randint(idx+1, len(nbrs)-1)
  return nbrs[i]

def higherDegNbrsCount(sortedAdjList, v):
  nbrMap, nbrs = sortedAdjList
  idx = nbrMap[v]
  return len(nbrs) - (idx + 1)

def sortedAdjLists(G):  
  sal = {}
  i = 0
  for u in G.Nodes():
    nbrs = [(G.GetNI(v).GetOutDeg(), v) for v in u.GetOutEdges()]
    nbrs = zip(*sorted(nbrs))[1]
    nbrMap = {nbr:i for i,nbr in enumerate(nbrs)}
    sal[u.GetId()] = (nbrMap, nbrs)
  return sal


def cumulativeEdgeDist2(G, sortedAdjLists):
  w = np.zeros(G.GetEdges())
  edgeMap = {}
  i = 0
  for e in G.Edges():
    u, v = e.GetId()
    if u == v:
      continue
    l_uv = higherDegNbrsCount(sortedAdjLists[u], v)
    l_vu = higherDegNbrsCount(sortedAdjLists[v], u)
    t = l_uv * l_vu
    w[i] = t
    edgeMap[i] = (u,v)
    i = i + 1
  W = sum(w)
  w1 = np.cumsum(w/W)
  return w1, W, edgeMap

def sample2(G, cumDist, edgeMap, sal):
  u, v = edgeMap[bisect.bisect(cumDist, random.random())]
  u1 = sampleHigherDegNbr(sal[u], v)
  v1 = sampleHigherDegNbr(sal[v], u)
  if u1 == v1: return None
  return (u1, u, v, v1)

def centeredSampler(G, k, cumDist, W, edgeMap, sal):
  #cumDist, W, edgeMap = cumulativeEdgeDist(G)
  mCount = [0]*6
  cBar = [0]*6
  i = 0
  for i in range(k):
    S_l = sample2(G, cumDist, edgeMap, sal)
    if S_l == None: continue
    u1, _, _, v1 = S_l
    if not G.IsEdge(u1, v1):
      continue
    mCount[inducedMotif(G, S_l)-1] += 1
    if i % 1000 == 0: print inducedMotif(G, S_l), i
  for i in range(4, 6+1):
    cBar[i-1] = (1.0 * mCount[i-1]/k) * W/A[4-1, i-1]
  #N_1 = sum([comb(n.GetOutDeg(), 3) for n in G.Nodes()])
  #cBar[1-1] = N_1 - cBar[3-1] - 2*cBar[5-1] - 4*cBar[6-1]
  return cBar


def compute4SubgraphFreq(G):
  w, W, edgeMap = cumulativeEdgeDist(G)
  nbrs = computeNbrs(G)
  cBar = threePathSampler(G, 200000, w, W, edgeMap, nbrs)

  sal = sortedAdjLists(G)
  w, W, edgeMap = cumulativeEdgeDist2(G, sal)
  cBar1 = centeredSampler(G, 200000, w, W, edgeMap, sal)
  for i in range(4, 6+1):
    cBar[i-1] = cBar1[i-1]

  return cBar


fileNames1 = ['higgs-mention_network.txt']
def computeSubgraphFreqs(topDir, subGraphSz, fileNames=None):
  start = os.getcwd()
  wf = open('aggr.freqs.' + str(subGraphSz), 'w+')
  for d, _, files in os.walk(topDir):
    os.chdir(d)
    print os.getcwd()
    gFiles = [f for f in files if f.split(".")[-1] == 'txt']
    for gFile in gFiles:
      if fileNames and gFile not in fileNames:
        continue
      G = snap.LoadEdgeList(snap.PUNGraph, gFile)
      if subGraphSz == 4:
        cBar = compute4SubgraphFreq(G)
        wf.write(gFile + ' ' +
                 ','.join([str(int(cBar[i])) for i in range(6)]) + '\n')
      else:
        ts,ps = countTriads(G)
        wf.write(gFile + ' ' + str(ps) + ',' + str(ts) + '\n')
        wf.flush()
        os.fsync(wf.fileno())
    os.chdir(topDir)
  os.chdir(start)
  wf.close()

computeSubgraphFreqs(os.getcwd() + '/data/rep/pa', 4)

