/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"

using namespace std;

// clustering

int network::hkFind(int x, vector<int>& labels) {
  int y = x;
  while (labels[y] != y) y = labels[y];

  while (labels[x] != x) {
    int z = labels[x];
    labels[x] = y;
    x = z;
  }
  return y;
}

int network::hkUnion(vector<int>& v, vector<int>& labels) {
  // joins equivalence classes and returns the canonical label of the resulting
  // class.
  for (unsigned i = 1; i < v.size(); ++i)
    labels[hkFind(v[i], labels)] = hkFind(v[0], labels);
  return hkFind(v[0], labels);
}

int network::hkMakeSet(vector<int>& labels) {
  // creates a new equivalence class and returns its label
  labels[0]++;
  labels.push_back(labels[0]);
  return labels[0];
}

void network::clusterPores(cluster* (pore::*getter)(void)const,
                           void (pore::*setter)(cluster*),
                           char (pore::*status)(void) const, char flag,
                           vector<cluster*>& clustersList) {
  if (!clustersList.empty())
    for (unsigned i = 0; i < clustersList.size(); ++i) delete clustersList[i];
  clustersList.clear();

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if ((p->*status)() == flag) p->setClusterTemp(0);
  }

  vector<int> labels;
  labels.push_back(0);

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if ((p->*status)() == flag) {
      vector<int> neighboorsClusters;
      vector<pore*> neighboors = p->getNeighboors();
      for (unsigned j = 0; j < neighboors.size(); j++) {
        if ((neighboors[j]->*status)() == flag &&
            neighboors[j]->getClusterTemp() != 0)
          neighboorsClusters.push_back(neighboors[j]->getClusterTemp());
      }
      if (neighboorsClusters.empty())
        p->setClusterTemp(hkMakeSet(labels));
      else if (neighboorsClusters.size() == 1)
        p->setClusterTemp(neighboorsClusters[0]);
      else
        p->setClusterTemp(hkUnion(neighboorsClusters, labels));
    }
  }
  /* apply the relabeling to the matrix */
  /* Create a mapping from the canonical labels determined by union/find into a
   * new set of canonical labels, which are guaranteed to be sequential. */

  vector<int> new_labels(labels.size(), 0);

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if ((p->*status)() == flag) {
      int x = hkFind(p->getClusterTemp(), labels);
      if (new_labels[x] == 0) {
        new_labels[0]++;
        new_labels[x] = new_labels[0];
        clustersList.push_back(new cluster(new_labels[0]));
      }
      (p->*setter)(clustersList[new_labels[x] - 1]);
    }
  }

  // Identify sepecial clusters
  set<cluster*> inletClusters, outletClusters, spanningClusters;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if ((p->*status)() == flag && p->getInlet())
      inletClusters.insert((p->*getter)());
    if ((p->*status)() == flag && p->getOutlet())
      outletClusters.insert((p->*getter)());
  }
  for (set<cluster*>::iterator iterator = inletClusters.begin();
       iterator != inletClusters.end(); ++iterator) {
    (*iterator)->setInlet(true);
  }
  for (set<cluster*>::iterator iterator = outletClusters.begin();
       iterator != outletClusters.end(); ++iterator) {
    (*iterator)->setOutlet(true);
  }
  std::set_intersection(
      inletClusters.begin(), inletClusters.end(), outletClusters.begin(),
      outletClusters.end(),
      std::inserter(spanningClusters, spanningClusters.begin()));
  for (set<cluster*>::iterator iterator = spanningClusters.begin();
       iterator != spanningClusters.end(); ++iterator) {
    (*iterator)->setSpanning(true);
  }
}

void network::clusterEverything() {
  cluster* (pore::*getter)() const = &pore::getClusterExist;
  void (pore::*setter)(cluster*) = &pore::setClusterExist;
  char (pore::*status)(void) const = &pore::getExist;
  clusterPores(getter, setter, status, 't', existClusters);
}
