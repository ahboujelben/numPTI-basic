/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef CLUSTER_H
#define CLUSTER_H

#include <map>
#include <vector>

class cluster {
 public:
  cluster(int);
  int getLabel() const;
  void setLabel(int);

  bool getInlet() const;
  void setInlet(bool value);

  bool getOutlet() const;
  void setOutlet(bool value);

  bool getSpanning() const;
  void setSpanning(bool value);

 private:
  int id;
  int label;
  bool inlet;
  bool outlet;
  bool spanning;
};

#endif  // CLUSTER_H
