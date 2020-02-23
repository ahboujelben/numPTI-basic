/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "element.h"

element::element() {
  radius = radius_sq = 0;
  length = length_sq = 0;
  viscosity = 1;
  volume = 0;
  conductivity = 0;
  concentration = 0;
  HDConcentration = 0;
  shapeFactor = 1. / (4 * 3.1415);
  shapeFactorConstant = 0.5;
  HDConcentration = 0;
  exist = 't';
  vesselType = 2;

  closed = false;
  inlet = false;
  outlet = false;

  clusterTemp = 0;
  clusterExist = 0;
}

element::~element() {}
int element::getId() const { return id; }

void element::setId(int value) { id = value; }
double element::getRadius() const { return radius; }

void element::setRadius(double value) { radius = value; }
double element::getLength() const { return length; }

void element::setLength(double value) { length = value; }
double element::getVolume() const { return volume; }

void element::setVolume(double value) { volume = value; }
double element::getShapeFactor() const { return shapeFactor; }

void element::setShapeFactor(double value) { shapeFactor = value; }
double element::getConductivity() const { return conductivity; }

void element::setConductivity(double value) { conductivity = value; }

int element::getClusterTemp() const { return clusterTemp; }

void element::setClusterTemp(int value) { clusterTemp = value; }

bool element::getClosed() const { return closed; }

void element::setClosed(bool value) {
  closed = value;
  exist = value == true ? 'f' : 't';
}
bool element::getInlet() const { return inlet; }

void element::setInlet(bool value) { inlet = value; }
bool element::getOutlet() const { return outlet; }

void element::setOutlet(bool value) { outlet = value; }

double element::getConcentration() const { return concentration; }

void element::setConcentration(double value) { concentration = value; }

double element::getViscosity() const { return viscosity; }

void element::setViscosity(double value) { viscosity = value; }
double element::getShapeFactorConstant() const { return shapeFactorConstant; }

void element::setShapeFactorConstant(double value) {
  shapeFactorConstant = value;
}

double element::getFlow() const { return flow; }

void element::setFlow(double value) { flow = value; }
double element::getMassFlow() const { return massFlow; }

void element::setMassFlow(double value) { massFlow = value; }

double element::getVesselType() const { return vesselType; }

void element::setVesselType(double value) { vesselType = value; }

cluster *element::getClusterExist() const { return clusterExist; }

void element::setClusterExist(cluster *value) { clusterExist = value; }
char element::getExist() const { return exist; }

void element::setExist(char value) { exist = value; }
double element::getRadius_sq() const { return radius_sq; }

void element::setRadius_sq(double value) { radius_sq = value; }
double element::getLength_sq() const { return length_sq; }

void element::setLength_sq(double value) { length_sq = value; }

double element::getHDConcentration() const { return HDConcentration; }

void element::setHDConcentration(double value) { HDConcentration = value; }

double element::getMembranePermeability() const { return membranePermeability; }

void element::setMembranePermeability(double value) {
  membranePermeability = value;
}
