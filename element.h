/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef ELEMENT_H
#define ELEMENT_H

#include "cluster.h"

class element {
 public:
  element();
  virtual ~element();

  int getId() const;
  void setId(int value);

  double getRadius() const;
  void setRadius(double value);

  double getLength() const;
  void setLength(double value);

  double getVolume() const;
  void setVolume(double value);

  double getShapeFactor() const;
  void setShapeFactor(double value);

  double getConductivity() const;
  void setConductivity(double value);

  int getClusterTemp() const;
  void setClusterTemp(int value);

  bool getClosed() const;
  void setClosed(bool value);

  bool getInlet() const;
  void setInlet(bool value);

  bool getOutlet() const;
  void setOutlet(bool value);

  double getConcentration() const;
  void setConcentration(double value);

  double getViscosity() const;
  void setViscosity(double value);

  double getShapeFactorConstant() const;
  void setShapeFactorConstant(double value);

  double getFlow() const;
  void setFlow(double value);

  double getMassFlow() const;
  void setMassFlow(double value);

  double getVesselType() const;
  void setVesselType(double value);

  char getExist() const;
  void setExist(char value);

  double getRadius_sq() const;
  void setRadius_sq(double value);

  double getLength_sq() const;
  void setLength_sq(double value);

  double getHDConcentration() const;
  void setHDConcentration(double value);

  double getMembranePermeability() const;
  void setMembranePermeability(double value);

  cluster *getClusterExist() const;
  void setClusterExist(cluster *value);

 protected:
  int type;

  int id;
  double radius;
  double radius_sq;
  double length;
  double length_sq;
  double volume;
  double shapeFactor;
  double shapeFactorConstant;
  double conductivity;
  double viscosity;
  double vesselType;
  double concentration;
  double HDConcentration;
  double flow;
  double massFlow;
  double membranePermeability;

  bool inlet;
  bool outlet;
  bool closed;
  char exist;

  int clusterTemp;
  cluster *clusterExist;
};

#endif  // ELEMENT_H
