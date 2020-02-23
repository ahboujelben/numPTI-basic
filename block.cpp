/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "block.h"

block::block(double X, double Y, double Z) {
  x = X;
  y = Y;
  z = Z;
  inlet = false;
  outlet = false;
  inletY = false;
  outletY = false;
  inletZ = false;
  outletZ = false;
  closed = false;
  volume = 0;
  conductivity = 0;
  concentration = 0;
  HDConcentration = 0;
  connectedToVessel = false;

  oxygenConcentration = 0;
  TAFConcentration = 0;
  MDEConcentration = 0;
  FNConcentration = 0;
}
int block::getId() const { return id; }

void block::setId(int value) { id = value; }
int block::getX() const { return x; }

void block::setX(int value) { x = value; }
int block::getY() const { return y; }

void block::setY(int value) { y = value; }
int block::getZ() const { return z; }

void block::setZ(int value) { z = value; }
double block::getXCoordinate() const { return xCoordinate; }

void block::setXCoordinate(double value) { xCoordinate = value; }
double block::getYCoordinate() const { return yCoordinate; }

void block::setYCoordinate(double value) { yCoordinate = value; }
double block::getZCoordinate() const { return zCoordinate; }

void block::setZCoordinate(double value) { zCoordinate = value; }
double block::getHx() const { return hx; }

void block::setHx(double value) { hx = value; }
double block::getHy() const { return hy; }

void block::setHy(double value) { hy = value; }
double block::getHz() const { return hz; }

void block::setHz(double value) { hz = value; }
double block::getVolume() const { return volume; }

void block::setVolume(double value) { volume = value; }
double block::getConductivity() const { return conductivity; }

void block::setConductivity(double value) { conductivity = value; }
double block::getConcentration() const { return concentration; }

void block::setConcentration(double value) { concentration = value; }
bool block::getInlet() const { return inlet; }

void block::setInlet(bool value) { inlet = value; }
bool block::getOutlet() const { return outlet; }

void block::setOutlet(bool value) { outlet = value; }
bool block::getInletY() const { return inletY; }

void block::setInletY(bool value) { inletY = value; }
bool block::getOutletY() const { return outletY; }

void block::setOutletY(bool value) { outletY = value; }
bool block::getInletZ() const { return inletZ; }

void block::setInletZ(bool value) { inletZ = value; }
bool block::getOutletZ() const { return outletZ; }

void block::setOutletZ(bool value) { outletZ = value; }
bool block::getClosed() const { return closed; }

void block::setClosed(bool value) { closed = value; }
std::vector<block *> block::getNeighbooringBlocks() const {
  return neighbooringBlocks;
}

void block::setNeighbooringBlocks(const std::vector<block *> &value) {
  neighbooringBlocks = value;
}

bool block::getConnectedToVessel() const { return connectedToVessel; }

void block::setConnectedToVessel(bool value) { connectedToVessel = value; }
double block::getEffectiveVolume() const { return effectiveVolume; }

void block::setEffectiveVolume(double value) { effectiveVolume = value; }
std::map<int, double> &block::neighbooringVesselsArea() {
  return neighbooringVesselsAreaMap;
}

double block::getOxygenConcentration() const { return oxygenConcentration; }

void block::setOxygenConcentration(double value) {
  oxygenConcentration = value;
}

double block::getTAFConcentration() const { return TAFConcentration; }

void block::setTAFConcentration(double value) { TAFConcentration = value; }

double block::getMDEConcentration() const { return MDEConcentration; }

void block::setMDEConcentration(double value) { MDEConcentration = value; }

double block::getFNConcentration() const { return FNConcentration; }

void block::setFNConcentration(double value) { FNConcentration = value; }
double block::getDiffusivity() const { return diffusivity; }

void block::setDiffusivity(double value) { diffusivity = value; }

double block::getHDConcentration() const { return HDConcentration; }

void block::setHDConcentration(double value) { HDConcentration = value; }
