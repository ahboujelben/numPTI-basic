/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef BLOCK_H
#define BLOCK_H

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

class block {
 public:
  block(double, double, double);

  int getId() const;
  void setId(int value);

  int getX() const;
  void setX(int value);

  int getY() const;
  void setY(int value);

  int getZ() const;
  void setZ(int value);

  double getXCoordinate() const;
  void setXCoordinate(double value);

  double getYCoordinate() const;
  void setYCoordinate(double value);

  double getZCoordinate() const;
  void setZCoordinate(double value);

  double getHx() const;
  void setHx(double value);

  double getHy() const;
  void setHy(double value);

  double getHz() const;
  void setHz(double value);

  double getVolume() const;
  void setVolume(double value);

  double getConductivity() const;
  void setConductivity(double value);

  double getConcentration() const;
  void setConcentration(double value);

  bool getInlet() const;
  void setInlet(bool value);

  bool getOutlet() const;
  void setOutlet(bool value);

  bool getInletY() const;
  void setInletY(bool value);

  bool getOutletY() const;
  void setOutletY(bool value);

  bool getInletZ() const;
  void setInletZ(bool value);

  bool getOutletZ() const;
  void setOutletZ(bool value);

  bool getClosed() const;
  void setClosed(bool value);

  bool getConnectedToVessel() const;
  void setConnectedToVessel(bool value);

  double getEffectiveVolume() const;
  void setEffectiveVolume(double value);

  double getOxygenConcentration() const;
  void setOxygenConcentration(double value);

  double getTAFConcentration() const;
  void setTAFConcentration(double value);

  double getMDEConcentration() const;
  void setMDEConcentration(double value);

  double getFNConcentration() const;
  void setFNConcentration(double value);

  double getDiffusivity() const;
  void setDiffusivity(double value);

  bool getColored() const;
  void setColored(bool value);

  double getHDConcentration() const;
  void setHDConcentration(double value);

  std::vector<block *> getNeighbooringBlocks() const;
  void setNeighbooringBlocks(const std::vector<block *> &value);

  std::map<int, double> &neighbooringVesselsArea();

 private:
  int id;
  int x;
  int y;
  int z;
  double xCoordinate;
  double yCoordinate;
  double zCoordinate;
  double hx;
  double hy;
  double hz;
  double volume;
  double conductivity;
  double concentration;
  double HDConcentration;
  double oxygenConcentration;
  double TAFConcentration;
  double MDEConcentration;
  double FNConcentration;
  double effectiveVolume;
  double diffusivity;

  bool inlet;
  bool outlet;
  bool inletY;
  bool outletY;
  bool inletZ;
  bool outletZ;
  bool connectedToVessel;

  bool closed;

  std::vector<block *> neighbooringBlocks;
  std::map<int, double> neighbooringVesselsAreaMap;
};

#endif  // BLOCK_H
