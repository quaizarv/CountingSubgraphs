import os
import collections
import numpy as np
import snap
import re
from scipy.misc import comb

#(degrees, clustering, average path length, size of giant component, degree corellation, assortativity

def ndDist(G):
  d = collections.Counter(n.GetOutDeg() for n in G.Nodes())
  Z = sum(d.values())
  ndd = collections.defaultdict(lambda: 0)
  for k,v in d.items():
    ndd[k] = v*1.0/Z
  return ndd

def avgDegreeOfNbrs(G, n):
  d = 0
  for nbr in n.GetOutEdges():
   d += G.GetNI(nbr).GetOutDeg()
  return d * 1.0/n.GetOutDeg()

def knnDist(G):
  knnTemp = collections.defaultdict(lambda: (0, 0))
  for n in G.Nodes():
    cnt, totalAvgDeg = knnTemp[n.GetOutDeg()]
    knnTemp[n.GetOutDeg()] = (cnt+1, totalAvgDeg + avgDegreeOfNbrs(G, n))

  knn = collections.defaultdict(lambda: 0)
  for k, (cnt, totalAvgDeg) in knnTemp.items():
    knn[k] = totalAvgDeg * 1.0/cnt

  return knn

def ccDist(G):
  ccTemp = collections.defaultdict(lambda: (0, 0))
  for n in G.Nodes():
    cnt, totalCC = ccTemp[n.GetOutDeg()]
    ccTemp[n.GetOutDeg()] = (cnt+1, totalCC + snap.GetNodeClustCf(G, n.GetId()))

  ccD = collections.defaultdict(lambda: 0)
  for k, (cnt, totalCC) in ccTemp.items():
    ccD[k] = totalCC * 1.0/cnt

  return ccD

def dkDist(G):
  dkD = collections.Counter()
  for e in G.Edges():
    u,v = e.GetId()
    dkD[(G.GetNI(u).GetOutDeg(), G.GetNI(v).GetOutDeg())] += 1
    dkD[(G.GetNI(v).GetOutDeg(), G.GetNI(u).GetOutDeg())] += 1
  return dkD

def ASCoeff(G):
  M = G.GetEdges()
  dudv, duPlusdv, du2Plusdv2 = (0, 0, 0)
  for e in G.Edges():
    u,v = e.GetId()
    du, dv = (G.GetNI(u).GetOutDeg(), G.GetNI(v).GetOutDeg())
    dudv += du * dv * 1.0/M
    duPlusdv += (du + dv) * 0.5/M
    du2Plusdv2 += (du**2 + dv**2) * 0.5/M
    
  return (dudv - duPlusdv**2)/(du2Plusdv2 - duPlusdv**2)

def euclidDist(v1, v2):
  dist = 0
  for i in set(v1.keys() + v2.keys()):
    dist += (v1[i] - v2[i])**2
  return dist ** 0.5

def computeGraphMetrics(GRep, GModel):
  nddDiff = euclidDist(ndDist(GRep), ndDist(GModel))
  knnDiff = euclidDist(knnDist(GRep), knnDist(GModel))
  dkDiff = euclidDist(dkDist(GRep), dkDist(GModel))
  ccDiff = euclidDist(ccDist(GRep), ccDist(GModel))
  ASVals = (ASCoeff(GRep), ASCoeff(GModel))
  MxWccVals = (snap.GetMxWcc(GRep).GetNodes(), snap.GetMxWcc(GModel).GetNodes())
  effDVals = (snap.GetBfsEffDiam(GRep, 1000, False), 
              snap.GetBfsEffDiam(GModel, 1000, False))
  return nddDiff, knnDiff, dkDiff, ccDiff, ASVals, MxWccVals, effDVals


def computeMetricsAllGraphs(topDir, fileNames=None):
  start = os.getcwd()
  for d, _, files in os.walk(topDir):
    os.chdir(d)
    print os.getcwd()
    gFiles = [f for f in files if f.split(".")[-1] == 'txt']
    print gFiles
    for gFile in gFiles:
      if fileNames and gFile not in fileNames:
        continue
      origGFile = '../' + gFile.split('.')[0] + '.txt'
      Gr = snap.LoadEdgeList(snap.PUNGraph, origGFile)
      Gm = snap.LoadEdgeList(snap.PUNGraph, gFile)
      print origGFile, gFile, Gr.GetNodes(), Gm.GetNodes()
      print computeGraphMetrics(Gr, Gm)      
    os.chdir(topDir)
  os.chdir(start)



origGFreqs4 = \
    {'ca-HepPh.txt': (143447956,204516498,462525447,886898,34986029,150346548),
     'twitter_combined.txt': (56145695401,16745491387,5749234503,138156186,571648938,105481650),
     'wiki-Vote.txt': (1127969256,1050479841,283523147,23135248,28104277,2072670),
     'loc-gowalla_edges.txt': (784255583714,15223729099,3114855610,41903138,85995370,6064492),
     'p2p-Gnutella31.txt': (8100733,15436049,69627,41717,750,14)}

nnGFreqs4 = \
    {'ca-HepPh.txt': (175398209,374023492,41337667,3321536,2417126,174761),
     'loc-gowalla_edges.txt': (352387326,674424500,46336606,610250,1273082,100746),
     'p2p-Gnutella31.txt': (95943568,216131490,12813389,300165,306610,15432),
     'twitter_combined.txt': (23589352203,38853134377,3962226150,277027959,210999159,15410609),
     'wiki-Vote.txt': (523299652,905388011,234517958,18187016,26383895,3306396)}

kronGFreqs4 = \
    {'ca-HepPh.txt': (301319862,129298742,6511985,389742,93707,1067),
     'loc-gowalla_edges.txt': (68483990501,19249636288,303416392,14304319,1331386,3318),
     'p2p-Gnutella31.txt': (4030000713,1443296154,36082403,1994202,247201,1246),
     'twitter_combined.txt': (19202934842,5904743497,132211714,6521889,851767,4440),
     'wiki-Vote.txt': (52731836,29336530,1416071,108473,20266,208)}

ffGFreqs4 = \
    {'ca-HepPh.txt': (70230227,55136139,27939585,981706,3690382,540356),
     'loc-gowalla_edges.txt': (256676088759,215909513611,45500533063,2559920149,3206439676,243908791),
     'p2p-Gnutella31.txt': (25583262293,18761233280,5728862866,255669089,502544277,51523126),
     'twitter_combined.txt': (884606946692,458236206376,466327186574,7848408553,80067923424,24599891944),
     'wiki-Vote.txt': (1146838155,450840505,751887014,6825659,171125618,62340757)}

paGFreqs4 = \
    {'ca-HepPh.txt': (22101450,8035619,145088,5103,934,10),
     'loc-gowalla_edges.txt': (4665455212,1051760272,5537903,101887,13678,179),
     'p2p-Gnutella31.txt': (137417186,22938714,109929,1643,191,0),
     'twitter_combined.txt': (21219324538,9807776590,175205280,5717167,1404717,30123),
     'wiki-Vote.txt': (352973374,319072672,17604134,1127633,432850,16413)}


def freqVecDistance(v1, v2):
  return np.linalg.norm(np.array(list(v1), dtype=np.double) - 
                        np.array(list(v2), dtype=np.double))

def readGraphSizeStats(fileName):
  graphSize = {}
  for line in file(fileName):
    items = re.split('\s+', line.strip())
    items = [item for item in items if item != '']
    if len(items) != 3:
      continue
    graphName = items[0]
    nodes = int(items[1])
    edges = int(items[2])
    graphSize[graphName] = (nodes, edges)
  return graphSize

def convertFreqCountsToProbs(freqStats, graphSize, subGraphSz):
  fs = {}
  for graphName, freqVec in freqStats.items():
    c = comb(graphSize[graphName][0], subGraphSz)
    temp = [elem * 1.0/c for elem in freqVec]
    #firstCoord = 1.0 - sum(temp)
    fs[graphName] = temp #[firstCoord] + temp
  return fs

def freqVecDistanceForAll(fRepMap, fModelMap, subGraphSz):
  distMap = {}
  graphSize = readGraphSizeStats('./results/graphSize.txt')
  frMap = convertFreqCountsToProbs(fRepMap, graphSize, subGraphSz)
  fmMap = convertFreqCountsToProbs(fModelMap, graphSize, subGraphSz)
  for n, v in frMap.items():
    distMap[n] = freqVecDistance(v, fmMap[n])
    print n, ':', 
    print '%2e'%distMap[n]
  #return distMap
