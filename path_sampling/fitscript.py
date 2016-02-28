import collections
import math
from subprocess import call
import numpy as np
import random
import snap

#GraphNames = ['p2p-Gnutella31']
GraphNames = ['loc-gowalla_edges']

GraphNames1 = ['as20000102',
              'wiki-Vote',
              'ca-HepPh',
              'ca-CondMat',
              'loc-brightkite_edges',
              'facebook_combined',
              'p2p-Gnutella31',
              'roadNet-CA'
               ]

nnCmd = 'python ./nearestNbr/generator/fitModel_noPickling.py -nn %s ./data/rep/%s.txt ./jdd.py /Users/qv/cs224w/fit/ -coarse'

ffCmd = 'python /Users/qv/cs224w/gd/project/sw/nearestNbr/generator/fitModel_noPickling.py -ff %s ./data/rep/%s.txt ./jdd.py /Users/qv/cs224w/fit/'

def fitPA():
  for gname in GraphNames:
    G = snap.LoadEdgeList(snap.PUNGraph, 'data/rep/' + gname + '.txt')
    G1 = snap.GenPrefAttach(G.GetNodes(), 
                            int(math.ceil(1.0 * G.GetEdges()/G.GetNodes())))
    snap.SaveEdgeList(G1, '/Users/qv/cs224w/fit/' + gname + '-PA' + '.txt')

def runFit(model):
  for gname in GraphNames:
    if model == 'ff':
      cmdStr = ffCmd % (gname, gname)
    else:
      cmdStr = nnCmd % (gname, gname)
    print cmdStr.split(' ')
    call(cmdStr.split(' '))
    #call(['python', 'jdd.py', './data/rep/as20000102.txt', 'test'])

runFit('ff')
#runFit('nn')


def normalize(theta):
  Z = sum(sum(theta))
  return (1.0/Z) * theta


def kroneckerGraphGen(a, b, c, d, m, edgeCount):
  idxChoices = [(0,0), (0,1), (1,0), (1,1)]
  def randQuadrant(theta):
    return idxChoices[np.random.choice(4, p=theta)]

  noiseLevel = 0.1
  theta = np.array([[a, b], [c, d]])

  if (np.diag(theta) == np.zeros(theta.shape[0])).any():
    raise ValueError("Initial graph must have all self-edges")
  theta = normalize(theta)
  theta_i = [None] * m
  for i in range(m):
    x_i = random.uniform(-noiseLevel, noiseLevel)
    theta_i[i] = theta + \
        np.array([[-2*x_i*a/(a + d), x_i], [x_i, -2*x_i*d/(a+d)]])
    theta_i[i].shape = (theta.size,)

  G = snap.TUNGraph.New()
  for u in range(2**m):
    G.AddNode(u)

  for j in range(edgeCount):
    if j % 10000 == 0: print j
    while True:
      u,v = (0,0)
      for i in range(m):
        r,c = randQuadrant(theta_i[i])
        u += r * (2 ** (m - 1 - i))
        v += c * (2 ** (m - 1 - i))
      if u != v and not G.IsEdge(u, v): 
        G.AddEdge(u, v)
        break

  return G


def logLikelihoodEstimate(data, x_m):
  data = [x for x in data if x >= x_m]
  n = len(data)
  return 1 + n/(sum([math.log(x) for x in data]) - n * math.log(x_m))

def findPrefAttachmentAlpha(G, x_m):
  d = collections.Counter(n.GetOutDeg() for n in G.Nodes())
  return logLikelihoodEstimate(d, x_m)

GraphSzAndAlpha = {'as20000102' : (34546, 420921, 2), 
                   'ca-CondMat' : (12008, 118521, 8),
                   'ca-HepPh' : (23133, 93497, 5),
                   'loc-brightkite_edges' : (58228, 214078, 2),
                   'p2p-Gnutella31' : (62586, 147892, 11)
                   }

def computeNNModelParams():
  nnModelParams = {}
  for gname, (n, e, x_m) in GraphSzAndAlpha.items():
    G = snap.LoadEdgeList(snap.PUNGraph, './data/rep/' + gname + '.txt')
    n = G.GetNodes()
    e = G.GetEdges()
    a = findPrefAttachmentAlpha(G, x_m)
    alpha = (e - n) * 1.0/n
    k = (5 - a) * alpha * (a - 1) / (4 * (1 + alpha * ((a - 1)**2)))
    beta = alpha / k
    u = beta / (1 + beta)
    T = e * n / (u * k)
    nnModelParams[gname] = u, k, T
  return nnModelParams
