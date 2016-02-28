import sys
import collections
import snap

def jointDegDist(G):
  jddMap = collections.Counter()
  for e in G.Edges():
    u,v = e.GetId()
    if u == v:
      continue
    uDeg = G.GetNI(u).GetOutDeg()
    vDeg = G.GetNI(v).GetOutDeg()
    """if (uDeg > vDeg):
      jddMap[(vDeg, uDeg)] += 1
    else:
      jddMap[(uDeg, vDeg)] += 1"""
    jddMap[(uDeg, vDeg)] += 1
    jddMap[(vDeg, uDeg)] += 1
  return jddMap


if __name__ == "__main__":
  print sys.argv
  assert(len(sys.argv) == 3) 
  inFileName = sys.argv[1]
  outFileName = sys.argv[2]
  G = snap.LoadEdgeList(snap.PUNGraph, inFileName)
  jddMap = jointDegDist(G)
  sortedJDD = sorted(jddMap.items())
  wf = open(outFileName, 'w+')
  for (d1, d2), cnt in sortedJDD:
    wf.write(' '.join([str(d1), str(d2), str(cnt)]) + '\n')
  wf.close()
