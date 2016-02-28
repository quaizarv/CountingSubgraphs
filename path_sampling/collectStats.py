import os
import re
import snap
from scipy.misc import comb
def readfreqStats(fileName, subGraphSz):
  coords = 2 if subGraphSz == 3 else 6
  freqStats = {}
  for line in file(fileName):
    items = re.split('\s+', line.strip())
    items = [item for item in items if item != '']
    if len(items) != 2:
      continue
    graphName = items[0]
    freqVec = [int(intStr) for intStr in items[1].split(',')]
    if len(freqVec) != coords:
      continue
    freqStats[graphName] = freqVec
  return freqStats

def GraphSizeStats(topDir, fileNames=None):
  graphStats = {}
  start = os.getcwd()
  for d, _, files in os.walk(topDir):
    os.chdir(d)
    #print os.getcwd()
    files = [f for f in files if f.split('.')[-1] == 'txt']
    for f in files:
      print f
      if fileNames and f not in fileNames:
        continue
      G = snap.LoadEdgeList(snap.PUNGraph, f)
      graphStats[f] = (G.GetNodes(), G.GetEdges())
    os.chdir(topDir)
  os.chdir(start)
  return graphStats

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

def collectStatsFromFile(filePath, subGraphSz):
  freqStats = readfreqStats(filePath + '/aggr.freqs.' + str(subGraphSz),
                            subGraphSz)
  graphSize = readGraphSizeStats(filePath + '/graphSize.txt')
  return convertFreqCountsToProbs(freqStats, graphSize, subGraphSz)

#collectStatsFromFiles(os.getcwd() + "/results/", 4)







