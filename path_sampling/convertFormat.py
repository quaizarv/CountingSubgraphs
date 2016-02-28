import collections
import re
import snap

def convertToGtrieFormat(fileName, outFileName):
  edges = []
  for line in file(fileName):
    items = re.split("\s", line)
    if items[0] == '#': continue
    u = int(items[0])
    v = int(items[1])
    edges.append((u, v))
  endPoints = zip(*edges)
  nodes = collections.Counter(endPoints[0] + endPoints[1])
  m = max(nodes.keys()) + 1
  nodeMap = None
  nodeMap = {}
  if m > len(nodes):
    nodeMap = {k:i+1 for i, k in enumerate(nodes.keys())}
  else:
    nodeMap = {i:i for i in range(len(nodes))}
    nodeMap[0] = m
  if nodeMap:
    for i in range(len(edges)):
      u, v = edges[i]
      edges[i] = (nodeMap[u], nodeMap[v])
  f = open(outFileName, 'w')
  for u,v in edges:
    f.write(str(u) + " " + str(v) + "\n")
  f.close()




