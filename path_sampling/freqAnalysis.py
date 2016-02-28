import math
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import collections
import pylab
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from collectStats import *

def freqClusterAnalysis(subGraphSz, numClus):
  freqStats = collectStatsFromFile(os.getcwd() + "/results", subGraphSz)
  idxToName = {i:n for i, n in enumerate(freqStats)}
  freqArray = [None]*len(freqStats)
  for i in range(len(freqStats)):
    freqArray[i] = freqStats[idxToName[i]]
  freqMtx = np.array(freqArray, dtype = np.double)
  freqMtx = normalize(freqMtx, axis=1)

  # For high dimension, use PCA to reduce dimensionality
  if subGraphSz == 4:
    pca = PCA(n_components=3)
    pca.fit(freqMtx)
    freqMtx = pca.transform(freqMtx)

  km = KMeans(n_clusters=numClus)
  km.fit_predict(freqMtx)
  err = sum([np.linalg.norm(km.cluster_centers_[km.labels_[i]] - freqMtx[i]) ** 2
             for i in range(len(idxToName))])

  numClus = km.cluster_centers_.shape[0]
  clusters = [{idxToName[i] for i in range(len(idxToName)) 
               if km.labels_[i] == c}
              for c in range(numClus)]
  repNames = zip(*[min([(np.linalg.norm(km.cluster_centers_[c] - freqMtx[i]), 
                         (idxToName[i], i))
                        for i in range(len(idxToName)) if km.labels_[i] == c]) 
                   for c in range(numClus)])[1]

  sReps = [sorted([(np.linalg.norm(km.cluster_centers_[c] - freqMtx[i]), 
                    (idxToName[i], i))
                   for i in range(len(idxToName)) if km.labels_[i] == c]) 
           for c in range(numClus)]

  repGs = ['p2p-Gnutella31.txt', 'wiki-Vote.txt', 'loc-gowalla_edges.txt',
           'ca-HepPh.txt', 'twitter_combined.txt']
  repNames = [(n, i) for i, n in idxToName.items() if n in repGs]
  if subGraphSz == 5:
    plotKmeansResult3D(freqMtx, km.labels_, repNames)
  else:
    plotKmeansResult2D(freqMtx, km.labels_, repNames)
  
  return clusters, repNames, err, sReps


"""plt.scatter(freqMtx[:, 0], freqMtx[:, 1], c=km.labels_)
_, label, (x, y) = tags[0]
plt.annotate(
            label, 
            xy = (x, y), xytext = (-20, 20),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
plt.show()
centers = [[0.02, 0.1], [0.05, 0.25], [0.2, 0.35], [0.6, 0.1]
km = KMeans(n_clusters=4, init = np.array(centers))
cr = zip(km.labels_, [idxToName[i] for i in range(len(idxToName))], freqMtx)
tags = [min([(np.linalg.norm(km.cluster_centers_[i] - m), n, m) for i, n,m in cr if i == c]) for c in range(8)]
err = sum([np.linalg.norm(km.cluster_centers_[i] - m) for i, n,m in cr])"""


def plotOutDegDistr(G):
  d = collections.Counter(n.GetOutDeg() for n in G.Nodes())
  plt.plot(d.keys(), [v * 1.0/sum(d.values()) for v in d.values()])
  plt.xscale('log')
  plt.yscale('log')
  plt.xlabel('out-degree')
  plt.ylabel('frequency')
  plt.show()


def logLikelihoodEstimate(data, x_m):
  data = [x for x in data if x >= x_m]
  n = len(data)
  return 1 + n/(sum([math.log(x) for x in data]) - n * math.log(x_m))

def findPrefAttachmentAlpha(G, x_m):
  d = collections.Counter(n.GetOutDeg() for n in G.Nodes())
  return logLikelihoodEstimate(d, x_m)


def plotKmeansResult2D(X, labels, reps):
  pylab.scatter(X[:, 1], X[:, 2], s = 100, c = labels)
  for c, (r, i) in enumerate(reps):
    _, x, y = X[i]
    if r in ('twitter_combined.txt', 'wiki-Vote.txt'):
      xyt = (-30, -30)
      vat = 'top'
    else:
      xyt = (-70, 20)
      vat = 'bottom'
    pylab.annotate(
      r.split('.')[0],
      xy = (x, y), xytext = xyt,
      textcoords = 'offset points', ha = 'right', va = vat,
      bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
      arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
  pylab.show()

def plotKmeansResult3D(X, labels, reps):
  """fig = plt.figure(1, figsize=(6, 4))
  plt.clf()
  ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
  plt.cla()
  ax.set_xlabel('0')
  ax.set_ylabel('1')
  ax.set_zlabel('2')
  ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=labels)"""

  fig = pylab.figure(figsize=(12, 9))
  ax = fig.add_subplot(111, projection = '3d')
  ax = Axes3D(fig, rect=[0, 0, .95, 1], elev=48, azim=134)
  ax.set_xlabel('1')
  ax.set_ylabel('0')
  ax.set_zlabel('2')
  ax.scatter(X[:, 1], X[:, 0], X[:, 2], c=labels, s = 100)

  for c, (r, i) in enumerate(reps):
    x, y, z = X[i]
    x2, y2, _ = proj3d.proj_transform(y,x,z, ax.get_proj())
    if r in ('twitter_combined.txt', 'wiki-Vote.txt'):
      xyt = (30, -30)
      vat = 'top'
    else:
      xyt = (-30, 30)
      vat = 'bottom'
    pylab.annotate(
      r.split('.')[0],
      xy = (x2, y2), xytext = xyt,
      textcoords = 'offset points', ha = 'right', va = vat,
      bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
      arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
  pylab.show()

