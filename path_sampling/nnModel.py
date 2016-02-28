import collections
import random
import snap
import math
import bisect

def genNNGraph(n, u, k, T):
  G = snap.TUNGraph.New()
  nodeCounter = 0 # Keeps track of ID of next available node

  degreeArray = [0 for i in range(0,n)]
  d = []
  N = [0] * 2

  while nodeCounter < n: #for j in range(T):
    if nodeCounter >= n: 
      break

    if random.random() < u:
      if len(d) == 0:
        continue

      x = random.choice(d) # Pick a node from list
      
      # Pick 2 unique nodes in the list of neighbors
      N = random.sample([nbr for nbr in G.GetNI(x).GetOutEdges()], 2) 

      if not G.IsEdge(N[0], N[1]): # If no edge exists between the 2, connect
        G.AddEdge(N[0], N[1])
        degreeArray[N[0]] += 1
        degreeArray[N[1]] += 1
        if degreeArray[N[0]] == 2:
          bisect.insort(d, N[0])

        if degreeArray[N[1]] == 2:
          bisect.insort(d, N[1])
    else:
      nodesSoFar = nodeCounter
      G.AddNode(nodeCounter)
      nodeCounter += 1
      if nodeCounter % 10000  == 0: print nodeCounter

      if nodeCounter == 1: # No use continuing if only one node in the graph
        continue

      a = random.randrange(0, nodesSoFar) # Pick a node in the graph
      G.AddEdge(a, nodesSoFar)
      degreeArray[a] += 1
      degreeArray[nodesSoFar] += 1

      if degreeArray[a] == 2:
        bisect.insort(d, a)
      if degreeArray[nodesSoFar] == 2:
        bisect.insort(d, nodesSoFar)

      for i in range(0, k): # Connect k random pairs in the graph
        N[0] = random.randint(0, nodeCounter - 1)
        N[1] = N[0]

        while N[1] == N[0]: # Ensure the two nodes are different (no self-loops)
          N[1] = random.randint(0, nodeCounter - 1)

          if not G.IsEdge(N[0], N[1]):
            G.AddEdge(N[0], N[1])

            degreeArray[N[0]] += 1
            degreeArray[N[1]] += 1

            if degreeArray[N[0]] == 2:
              bisect.insort(d, N[0])
            if degreeArray[N[1]] == 2:
              bisect.insort(d, N[1])

  return G


def logLikelihoodEstimate(data, x_m):
  data = [x for x in data if x >= x_m]
  n = len(data)
  return 1 + n/(sum([math.log(x) for x in data]) - n * math.log(x_m))

def findPrefAttachmentAlpha(G, x_m):
  d = collections.Counter(n.GetOutDeg() for n in G.Nodes())
  return logLikelihoodEstimate(d, x_m)

def computeNNModelParams():
  nnModelParams = {}
  for gname, (_, _, x_m) in GraphSzAndAlpha.items():
    G = snap.LoadEdgeList(snap.PUNGraph, './data/rep/' + gname + '.txt')
    n = G.GetNodes()
    e = G.GetEdges()
    gamma = findPrefAttachmentAlpha(G, x_m)
    print gname, n, e, gamma
    beta = 1.0/(gamma - 1.0)
    a = 4.0*beta*beta - 4.0*beta + 4.0*e/n + 4.0
    b =  4.0*beta - 8.0*beta*beta - 4.0*e/n
    c = 4.0*beta*beta
    print b**2, 4*a*c
    u1 = (-b + math.sqrt(b*b - 4*a*c))/(2*a)
    u2 = (-b - math.sqrt(b*b - 4*a*c))/(2*a)
    assert(u1 > 0 or u2 > 0)
    assert(u1 < 1 or u2 < 1)
    if u1 > 0 and u1 < 1 and (u2 > 1 or u2 < 0): u = u1
    elif u2 > 0 and u2 < 1 and (u1 > 1 or u1 < 0): u = u2
    else: u = max(u1, u2)
    k = (1.0/n) * (e - n/(1-u))
    T = n /(1-u)
    print u, k, T
    beta = (u/(2*(1-u))) * ( -1 + math.sqrt(1 + 4*(k + 1)*(1-u)/u))
    print (1 - u)*T, n*k +  T, 1 + 1.0/beta
    nnModelParams[gname] = u, k, T
  return nnModelParams

GraphSzAndAlpha = {#'as20000102' : (6474, 13895, 2), 
                   #'ca-CondMat' : (23133, 93497, 10),
                   'wiki-Vote': (7115, 100762, 10),
                   'ca-HepPh' : (12008, 118521, 8),
                   #'loc-gowalla_edges': (196591, 950327, 20),
                   #'loc-brightkite_edges' : (58228, 214078, 20),
                   'twitter_combined': (81306, 1342310, 30),
                   #'p2p-Gnutella31' : (62586, 147892, 11)
                   }

