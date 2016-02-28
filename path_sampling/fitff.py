
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
  return sorted(jddMap.items)

def fitForestFire(G):
  for i in range(1, 11):
    p = 0.1 * i
    genG = snap.GenForestFire(G.GetNodes(), p, 0)
