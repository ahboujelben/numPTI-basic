/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "pore.h"
using namespace std;

pore::pore(node *const &pNodeIn, node *const &pNodeOut) {
  nodeIn = pNodeIn;
  nodeOut = pNodeOut;

  fullLength = 0;
  nodeInLength = 0;
  nodeOutLength = 0;
  parentVessel = false;
  visited = true;
}

pore::~pore() {}

std::vector<pore *> &pore::getNeighboors() { return neighboors; }

void pore::setNeighboors(const std::vector<pore *> &value) {
  neighboors = value;
}

double pore::getNodeInLength() const { return nodeInLength; }

void pore::setNodeInLength(double value) { nodeInLength = value; }
double pore::getNodeOutLength() const { return nodeOutLength; }

void pore::setNodeOutLength(double value) { nodeOutLength = value; }

double pore::getMinXCoordinate() const {
  if (nodeIn == 0) return nodeOut->getXCoordinate();
  if (nodeOut == 0) return nodeIn->getXCoordinate();
  return min(nodeIn->getXCoordinate(), nodeOut->getXCoordinate());
}

double pore::getMaxXCoordinate() const {
  if (nodeIn == 0) return nodeOut->getXCoordinate();
  if (nodeOut == 0) return nodeIn->getXCoordinate();
  return max(nodeIn->getXCoordinate(), nodeOut->getXCoordinate());
}

double pore::getMinYCoordinate() const {
  if (nodeIn == 0) return nodeOut->getYCoordinate();
  if (nodeOut == 0) return nodeIn->getYCoordinate();
  return min(nodeIn->getYCoordinate(), nodeOut->getYCoordinate());
}

double pore::getMaxYCoordinate() const {
  if (nodeIn == 0) return nodeOut->getYCoordinate();
  if (nodeOut == 0) return nodeIn->getYCoordinate();
  return max(nodeIn->getYCoordinate(), nodeOut->getYCoordinate());
}

double pore::getMinZCoordinate() const {
  if (nodeIn == 0) return nodeOut->getZCoordinate();
  if (nodeOut == 0) return nodeIn->getZCoordinate();
  return min(nodeIn->getZCoordinate(), nodeOut->getZCoordinate());
}

double pore::getMaxZCoordinate() const {
  if (nodeIn == 0) return nodeOut->getZCoordinate();
  if (nodeOut == 0) return nodeIn->getZCoordinate();
  return max(nodeIn->getZCoordinate(), nodeOut->getZCoordinate());
}

double pore::getXCoordinate() const {
  return (getMaxXCoordinate() + getMinXCoordinate()) / 2;
}

double pore::getYCoordinate() const {
  return (getMaxYCoordinate() + getMinYCoordinate()) / 2;
}

double pore::getZCoordinate() const {
  return (getMaxZCoordinate() + getMinZCoordinate()) / 2;
}

double pore::getFullLength() const { return fullLength; }

void pore::setFullLength(double value) { fullLength = value; }

std::map<int, double> &pore::neighbooringBlocksArea() {
  return neighbooringBlocksAreaMap;
}

void pore::setTissueNeighboorsArea(const std::map<int, double> &value) {
  neighbooringBlocksAreaMap = value;
}

std::map<int, double> &pore::getFeedingVessels() { return feedingVessels; }

double pore::getInflowShare() const { return inflowShare; }

void pore::setInflowShare(double value) { inflowShare = value; }

bool pore::getVisited() const { return visited; }

void pore::setVisited(bool value) { visited = value; }
double pore::getAveragePressure() const { return averagePressure; }

void pore::setAveragePressure(double value) { averagePressure = value; }
double pore::getWSS() const { return WSS; }

void pore::setWSS(double value) { WSS = value; }
bool pore::getParentVessel() const { return parentVessel; }

void pore::setParentVessel(bool value) { parentVessel = value; }
double pore::getFQE() const { return FQE; }

void pore::setFQE(double value) { FQE = value; }

double pore::getConvectedStim() const { return convectedStim; }

void pore::setConvectedStim(double value) { convectedStim = value; }

double pore::getConductedStim() const { return conductedStim; }

void pore::setConductedStim(double value) { conductedStim = value; }

node *pore::getNodeIn() const { return nodeIn; }

void pore::setNodeIn(node *value) { nodeIn = value; }
node *pore::getNodeOut() const { return nodeOut; }

void pore::setNodeOut(node *value) { nodeOut = value; }
