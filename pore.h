/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef PORE_H
#define PORE_H

#include "cluster.h"
#include "element.h"
#include "node.h"

#include <cmath>
#include <vector>

class pore : public element {
 public:
  pore(node *const &, node *const &);
  ~pore();

  node *getNodeIn() const;
  void setNodeIn(node *value);
  node *getNodeOut() const;
  void setNodeOut(node *value);

  double getNodeInLength() const;
  void setNodeInLength(double value);

  double getNodeOutLength() const;
  void setNodeOutLength(double value);

  double getMinXCoordinate() const;
  double getMaxXCoordinate() const;

  double getMinYCoordinate() const;
  double getMaxYCoordinate() const;

  double getMinZCoordinate() const;
  double getMaxZCoordinate() const;

  double getXCoordinate() const;
  double getYCoordinate() const;
  double getZCoordinate() const;

  double getFullLength() const;
  void setFullLength(double value);

  double getInflowShare() const;
  void setInflowShare(double value);

  bool getVisited() const;
  void setVisited(bool value);

  double getAveragePressure() const;
  void setAveragePressure(double value);

  double getWSS() const;
  void setWSS(double value);

  bool getParentVessel() const;
  void setParentVessel(bool value);

  double getFQE() const;
  void setFQE(double value);

  double getConvectedStim() const;
  void setConvectedStim(double value);

  double getConductedStim() const;
  void setConductedStim(double value);

  std::vector<pore *> &getNeighboors();
  void setNeighboors(const std::vector<pore *> &value);

  std::map<int, double> &neighbooringBlocksArea();
  void setTissueNeighboorsArea(const std::map<int, double> &value);

  std::map<int, double> &getFeedingVessels();

 protected:
  node *nodeIn;
  node *nodeOut;

  double fullLength;
  double nodeInLength;
  double nodeOutLength;

  std::vector<pore *> neighboors;

  // Drug flow attributes
  std::map<int, double> feedingVessels;
  std::map<int, double> neighbooringBlocksAreaMap;

  // Angiogenesis attributes
  double WSS;
  double FQE;
  double convectedStim;
  double conductedStim;
  double inflowShare;
  double averagePressure;
  bool parentVessel;
  bool visited;
};

#endif  // PORE_H
