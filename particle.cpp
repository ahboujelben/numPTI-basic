/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
////////////////////////////////////////////////////////////////////////////

#include "particle.h"

particle::particle() {
  id = 0;
  poreID = 0;
  xCoordinate = 0;
  yCoordinate = 0;
  zCoordinate = 0;
  porePosition = 0;
  closed = false;
}
int particle::getId() const { return id; }

void particle::setId(int value) { id = value; }
int particle::getPoreID() const { return poreID; }

void particle::setPoreID(int value) { poreID = value; }
double particle::getXCoordinate() const { return xCoordinate; }

void particle::setXCoordinate(double value) { xCoordinate = value; }
double particle::getYCoordinate() const { return yCoordinate; }

void particle::setYCoordinate(double value) { yCoordinate = value; }
double particle::getZCoordinate() const { return zCoordinate; }

void particle::setZCoordinate(double value) { zCoordinate = value; }
double particle::getPorePosition() const { return porePosition; }

void particle::setPorePosition(double value) { porePosition = value; }
bool particle::getClosed() const { return closed; }

void particle::setClosed(bool value) { closed = value; }
