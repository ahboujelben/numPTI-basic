/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "cluster.h"

cluster::cluster(int pLabel) {
  label = pLabel;
  inlet = false;
  outlet = false;
  spanning = false;
}
int cluster::getLabel() const { return label; }
void cluster::setLabel(int pLabel) { label = pLabel; }
bool cluster::getInlet() const { return inlet; }

void cluster::setInlet(bool value) { inlet = value; }
bool cluster::getOutlet() const { return outlet; }

void cluster::setOutlet(bool value) { outlet = value; }
bool cluster::getSpanning() const { return spanning; }

void cluster::setSpanning(bool value) { spanning = value; }
