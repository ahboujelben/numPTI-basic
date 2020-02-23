/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"

using namespace std;

double network::runHaematocritFlow() {
  double hx = xEdgeLength / meshSizeX;
  double hy = yEdgeLength / meshSizeY;
  double hz = zEdgeLength / meshSizeZ;

  double coefX = 1 / pow(hx, 2);
  double coefY = 1 / pow(hy, 2);
  double coefZ = 1 / pow(hz, 2);
  double coeff = 1 / pow(hx, 2) + 1 / pow(hy, 2) + 1 / pow(hz, 2);

  // set time step

  double flowTimeStep = 1e50;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && abs(p->getFlow()) > 1e-20) {
      double sumSource = 0;
      auto& neigh = p->neighbooringBlocksArea();
      for (auto iterator : neigh) {
        int blockID = iterator.first;
        block* bb = getBlock(blockID - 1);
        if (!bb->getClosed()) {
          double area = iterator.second;
          sumSource += p->getMembranePermeability() * area / (p->getVolume());
        }
      }

      double step = 1. / (abs(p->getFlow()) / p->getVolume() + sumSource);
      if (step < flowTimeStep) flowTimeStep = step;
    }
  }

  if (updateBlockAttributes)
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        double sumSource = 0;

        map<int, double>& neigh = n->neighbooringVesselsArea();
        for (auto iterator : neigh) {
          int poreID = iterator.first;
          pore* pp = getPore(poreID - 1);
          if (!pp->getClosed()) {
            double area = iterator.second;
            sumSource +=
                pp->getMembranePermeability() * area / (n->getVolume());
          }
        }

        double step =
            1. / (2 * n->getDiffusivity() * coeff + sumSource + sigma);
        if (step < flowTimeStep) flowTimeStep = step;
      }
    }

  double deltaT = flowTimeStep;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      if (!p->getInlet())  // inlet pore
      {
        p->getFeedingVessels().clear();
        double totalInflow(0);
        vector<pore*> connectedPores = p->getNeighboors();
        for (unsigned j = 0; j < connectedPores.size(); j++) {
          pore* npore = connectedPores[j];
          if (!npore->getClosed()) {
            if ((p->getFlow() > 0 && npore->getFlow() > 0 &&
                 npore->getNodeIn() == p->getNodeOut()) ||
                (p->getFlow() > 0 && npore->getFlow() < 0 &&
                 npore->getNodeOut() == p->getNodeOut()) ||
                (p->getFlow() < 0 && npore->getFlow() > 0 &&
                 npore->getNodeIn() == p->getNodeIn()) ||
                (p->getFlow() < 0 && npore->getFlow() < 0 &&
                 npore->getNodeOut() == p->getNodeIn())) {
              p->getFeedingVessels()[npore->getId()] = abs(npore->getFlow());
              totalInflow += abs(npore->getFlow());
            }
          }
        }
        p->setInflowShare(abs(p->getFlow()) / totalInflow);
      }
    }
  }

  vector<double> blockConcentration;
  vector<double> poreConcentration;

  if (updateBlockAttributes) blockConcentration.reserve(totalBlocks);
  poreConcentration.reserve(totalPores);

  if (updateBlockAttributes)
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        double ii = n->getX();
        double jj = n->getY();
        double kk = n->getZ();

        block *nW, *nE, *nN, *nS, *nU, *nD;
        nW = getBlock(ii - 1, jj, kk);
        nE = getBlock(ii + 1, jj, kk);
        nS = getBlock(ii, jj - 1, kk);
        nN = getBlock(ii, jj + 1, kk);
        nD = getBlock(ii, jj, kk - 1);
        nU = getBlock(ii, jj, kk + 1);

        double sumSource = 0;
        double sumSource2 = 0;

        map<int, double>& neigh = n->neighbooringVesselsArea();
        for (auto iterator : neigh) {
          int poreID = iterator.first;
          pore* pp = getPore(poreID - 1);
          if (!pp->getClosed()) {
            double area = iterator.second;
            sumSource +=
                pp->getMembranePermeability() * area / (n->getVolume());
            sumSource2 += pp->getHDConcentration() *
                          pp->getMembranePermeability() * area /
                          (n->getVolume());
          }
        }

        double newConcentration(0);

        newConcentration += n->getHDConcentration() *
                            (1 - deltaT * (2 * n->getDiffusivity() * coeff +
                                           sumSource + sigma));
        if (nW != 0)
          if (!nW->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * nW->getHDConcentration();
        if (nE != 0)
          if (!nE->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * nE->getHDConcentration();
        if (nS != 0)
          if (!nS->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * nS->getHDConcentration();
        if (nN != 0)
          if (!nN->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * nN->getHDConcentration();
        if (nD != 0)
          if (!nD->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * nD->getHDConcentration();
        if (nU != 0)
          if (!nU->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * nU->getHDConcentration();

        if (closedBoundaries) {
          if (nW == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * n->getHDConcentration();
          if (nE == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * n->getHDConcentration();
          if (nS == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * n->getHDConcentration();
          if (nN == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * n->getHDConcentration();
          if (nD == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * n->getHDConcentration();
          if (nU == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * n->getHDConcentration();
        }

        newConcentration += deltaT * sumSource2;

        blockConcentration.push_back(newConcentration);
      }
    }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double sumSource = 0;
      double sumSource2 = 0;
      double sumInflow = 0;

      auto& neigh = p->neighbooringBlocksArea();
      for (auto iterator : neigh) {
        int blockID = iterator.first;
        block* bb = getBlock(blockID - 1);
        if (!bb->getClosed()) {
          double area = iterator.second;
          sumSource += p->getMembranePermeability() * area / (p->getVolume());
          sumSource2 += bb->getHDConcentration() *
                        p->getMembranePermeability() * area / (p->getVolume());
        }
      }

      if (p->getInlet()) {
        double concentration = 1;
        sumInflow = concentration * abs(p->getFlow()) / p->getVolume();
      } else {
        auto& neighPores = p->getFeedingVessels();
        for (auto iterator : neighPores) {
          pore* pp = getPore(iterator.first - 1);
          double inflow = iterator.second;
          if (!pp->getClosed()) {
            sumInflow += pp->getHDConcentration() * p->getInflowShare() *
                         inflow / p->getVolume();
          }
        }
      }

      double newConcentration(0);
      newConcentration +=
          p->getHDConcentration() *
          (1 - deltaT * (abs(p->getFlow()) / p->getVolume() + sumSource));
      newConcentration += deltaT * sumInflow;
      newConcentration += deltaT * sumSource2;

      poreConcentration.push_back(newConcentration);
    }
  }

  unsigned j(0);
  if (updateBlockAttributes)
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        n->setHDConcentration(blockConcentration[j]);
        j++;
        if (n->getHDConcentration() < 0 || n->getHDConcentration() > 1) {
          cout << "block concentration out of range: "
               << n->getHDConcentration() << endl;
          cancel = true;
        }
      }
    }

  j = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setHDConcentration(poreConcentration[j]);
      j++;
      if (p->getHDConcentration() < 0 || p->getHDConcentration() > 1.001) {
        cout << "pore concentration out of range: " << p->getHDConcentration()
             << endl;
        cancel = true;
      }
    }
  }

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      double concentration(0);
      double neighboorsNumber(0);
      vector<int> connectedPores = n->getConnectedPores();
      for (unsigned j = 0; j < connectedPores.size(); j++) {
        pore* npore = getPore(connectedPores[j] - 1);
        if (!npore->getClosed()) {
          concentration += npore->getHDConcentration();
          neighboorsNumber++;
        }
      }
      n->setHDConcentration(concentration / neighboorsNumber);
    }
  }

  emitPlotSignal();

  return flowTimeStep * flowRate;
}

double network::runHaematocritFlowWithPhaseSeparation() {
  double hx = xEdgeLength / meshSizeX;
  double hy = yEdgeLength / meshSizeY;
  double hz = zEdgeLength / meshSizeZ;

  double coefX = 1 / pow(hx, 2);
  double coefY = 1 / pow(hy, 2);
  double coefZ = 1 / pow(hz, 2);
  double coeff = 1 / pow(hx, 2) + 1 / pow(hy, 2) + 1 / pow(hz, 2);

  // set time step

  double flowTimeStep = 1e50;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && abs(p->getFlow()) > 1e-40 &&
        p->getConductivity() != 1e-200) {
      double sumSource = 0;
      auto& neigh = p->neighbooringBlocksArea();
      for (auto iterator : neigh) {
        int blockID = iterator.first;
        block* bb = getBlock(blockID - 1);
        if (!bb->getClosed()) {
          double area = iterator.second;
          sumSource += p->getMembranePermeability() * area / (p->getVolume());
        }
      }

      double step = 1. / (abs(p->getFlow()) / p->getVolume() + sumSource);
      if (step < flowTimeStep) flowTimeStep = step;
    }
  }

  if (flowTimeStep == 1e50) return 0;

  if (updateBlockAttributes)
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        double sumSource = 0;

        map<int, double>& neigh = n->neighbooringVesselsArea();
        for (auto iterator : neigh) {
          int poreID = iterator.first;
          pore* pp = getPore(poreID - 1);
          if (!pp->getClosed()) {
            double area = iterator.second;
            sumSource +=
                pp->getMembranePermeability() * area / (n->getVolume());
          }
        }

        double step =
            1. / (2 * n->getDiffusivity() * coeff + sumSource + sigma);
        if (step < flowTimeStep) flowTimeStep = step;
      }
    }

  double deltaT = flowTimeStep;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      if (!p->getInlet())  // inlet pore
      {
        p->getFeedingVessels().clear();
        double totalInflow(0);
        vector<pore*> connectedPores = p->getNeighboors();
        for (unsigned j = 0; j < connectedPores.size(); j++) {
          pore* npore = connectedPores[j];
          if (!npore->getClosed()) {
            if ((p->getFlow() > 0 && npore->getFlow() > 0 &&
                 npore->getNodeIn() == p->getNodeOut()) ||
                (p->getFlow() > 0 && npore->getFlow() < 0 &&
                 npore->getNodeOut() == p->getNodeOut()) ||
                (p->getFlow() < 0 && npore->getFlow() > 0 &&
                 npore->getNodeIn() == p->getNodeIn()) ||
                (p->getFlow() < 0 && npore->getFlow() < 0 &&
                 npore->getNodeOut() == p->getNodeIn())) {
              p->getFeedingVessels()[npore->getId()] = abs(npore->getFlow());
              totalInflow += abs(npore->getFlow());
            }
          }
        }
        p->setInflowShare(abs(p->getFlow()) / totalInflow);
      }
    }
  }

  vector<double> blockConcentration;
  vector<double> poreConcentration;

  if (updateBlockAttributes) blockConcentration.reserve(totalBlocks);
  poreConcentration.reserve(totalPores);

  if (updateBlockAttributes)
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        double ii = n->getX();
        double jj = n->getY();
        double kk = n->getZ();

        block *nW, *nE, *nN, *nS, *nU, *nD;
        nW = getBlock(ii - 1, jj, kk);
        nE = getBlock(ii + 1, jj, kk);
        nS = getBlock(ii, jj - 1, kk);
        nN = getBlock(ii, jj + 1, kk);
        nD = getBlock(ii, jj, kk - 1);
        nU = getBlock(ii, jj, kk + 1);

        double sumSource = 0;
        double sumSource2 = 0;

        map<int, double>& neigh = n->neighbooringVesselsArea();
        for (auto iterator : neigh) {
          int poreID = iterator.first;
          pore* pp = getPore(poreID - 1);
          if (!pp->getClosed()) {
            double area = iterator.second;
            sumSource +=
                pp->getMembranePermeability() * area / (n->getVolume());
            sumSource2 += pp->getHDConcentration() *
                          pp->getMembranePermeability() * area /
                          (n->getVolume());
          }
        }

        double newConcentration(0);

        newConcentration += n->getHDConcentration() *
                            (1 - deltaT * (2 * n->getDiffusivity() * coeff +
                                           sumSource + sigma));
        if (nW != 0)
          if (!nW->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * nW->getHDConcentration();
        if (nE != 0)
          if (!nE->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * nE->getHDConcentration();
        if (nS != 0)
          if (!nS->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * nS->getHDConcentration();
        if (nN != 0)
          if (!nN->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * nN->getHDConcentration();
        if (nD != 0)
          if (!nD->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * nD->getHDConcentration();
        if (nU != 0)
          if (!nU->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * nU->getHDConcentration();

        if (closedBoundaries) {
          if (nW == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * n->getHDConcentration();
          if (nE == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * n->getHDConcentration();
          if (nS == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * n->getHDConcentration();
          if (nN == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * n->getHDConcentration();
          if (nD == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * n->getHDConcentration();
          if (nU == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * n->getHDConcentration();
        }

        newConcentration += deltaT * sumSource2;

        blockConcentration.push_back(newConcentration);
      }
    }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double sumSource = 0;
      double sumSource2 = 0;
      double sumInflow = 0;

      auto& neigh = p->neighbooringBlocksArea();
      for (auto iterator : neigh) {
        int blockID = iterator.first;
        block* bb = getBlock(blockID - 1);
        if (!bb->getClosed()) {
          double area = iterator.second;
          sumSource += p->getMembranePermeability() * area / (p->getVolume());
          sumSource2 += bb->getHDConcentration() *
                        p->getMembranePermeability() * area / (p->getVolume());
        }
      }

      if (p->getInlet()) {
        double concentration = 1;
        sumInflow = concentration * abs(p->getFlow()) / p->getVolume();
      } else {
        auto& neighPores = p->getFeedingVessels();
        for (auto iterator : neighPores) {
          pore* pp = getPore(iterator.first - 1);
          double inflow = iterator.second;
          if (!pp->getClosed()) {
            sumInflow += pp->getHDConcentration() * p->getFQE() * inflow /
                         p->getVolume();
          }
        }
      }

      double newConcentration(0);
      newConcentration +=
          p->getHDConcentration() *
          (1 - deltaT * (abs(p->getFlow()) / p->getVolume() + sumSource));
      newConcentration += deltaT * sumInflow;
      newConcentration += deltaT * sumSource2;

      poreConcentration.push_back(newConcentration);
    }
  }

  unsigned j(0);
  if (updateBlockAttributes)
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        n->setHDConcentration(blockConcentration[j]);
        j++;
        if (n->getHDConcentration() < 0 || n->getHDConcentration() > 1) {
          cout << "block concentration out of range: "
               << n->getHDConcentration() << endl;
          cancel = true;
        }
      }
    }

  j = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setHDConcentration(poreConcentration[j]);
      j++;
      if (poreConcentration[i] > 1) p->setHDConcentration(1);
      if (p->getHDConcentration() < 0 || p->getHDConcentration() > 1.001) {
        cout << "pore concentration out of range: " << p->getHDConcentration()
             << endl;
        cancel = true;
      }
    }
  }

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      double concentration(0);
      double neighboorsNumber(0);
      vector<int> connectedPores = n->getConnectedPores();
      for (unsigned j = 0; j < connectedPores.size(); j++) {
        pore* npore = getPore(connectedPores[j] - 1);
        if (!npore->getClosed()) {
          concentration += npore->getHDConcentration();
          neighboorsNumber++;
        }
      }
      n->setHDConcentration(concentration / neighboorsNumber);
    }
  }

  emitPlotSignal();

  return flowTimeStep * flowRate;
}

void network::assignInitialBloodViscosities() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double plasmaViscosity = 1.2e-3;
      double R = p->getRadius() * 1e6;
      double HD = p->getHDConcentration() * 0.45;
      double C =
          (0.8 + exp(-0.15 * R)) * (-1 + 1 / (1 + 1e-11 * pow(2 * R, 12))) +
          1 / (1 + 1e-11 * pow(2 * R, 12));
      double fHD = (pow(1 - HD, C) - 1) / (pow(1 - 0.45, C) - 1);
      double mu045 =
          6 * exp(-0.17 * R) + 3.2 - 2.44 * exp(-0.06 * pow(2 * R, 0.645));
      double murel = (1 + (mu045 - 1) * fHD * pow(2 * R / (2 * R - 1.1), 2)) *
                     pow(2 * R / (2 * R - 1.1), 2);
      double newViscosity = plasmaViscosity * murel;
      p->setViscosity(newViscosity);
    }
  }
}

void network::assignBloodViscosities() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && !p->getParentVessel()) {
      double R = p->getRadius();
      double HD = p->getHDConcentration() * 0.45;
      double C =
          (0.8 + exp(-0.15 * R)) * (-1 + 1 / (1 + 1e-11 * pow(2 * R, 12))) +
          1 / (1 + 1e-11 * pow(2 * R, 12));
      double fHD = (pow(1 - HD, C) - 1) / (pow(1 - 0.45, C) - 1);
      double mu045 = 6 * exp(-0.17 * R * 1e6) + 3.2 -
                     2.44 * exp(-0.06 * pow(2 * R * 1e6, 0.645));
      double murel =
          (1 + (mu045 - 1) * fHD * pow(2 * R * 1e6 / (2 * R * 1e6 - 1.1), 2)) *
          pow(2 * R * 1e6 / (2 * R * 1e6 - 1.1), 2);
      double newViscosity = plasmaViscosity * murel;
      p->setViscosity(newViscosity);
    }
  }
}

void network::calculateConvectedStimuli() {
  double nutrientRef = QHDref * 0.45;

  for (int j = 0; j < totalPores; ++j) {
    pore* p = getPore(j);
    if (!p->getClosed()) p->setConvectedStim(0);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && abs(p->getFlow()) > 1e-20) {
      double localNutrient = abs(p->getFlow()) * p->getHDConcentration() * 0.45;
      if (localNutrient < nutrientRef) {
        for (int j = 0; j < totalPores; ++j) {
          pore* p = getPore(j);
          if (!p->getClosed()) p->setVisited(false);
        }

        double addedStimulus =
            p->getLength() * 1.0e6 * (1 - localNutrient / nutrientRef);
        p->setConvectedStim(p->getConvectedStim() + addedStimulus);
        p->setVisited(true);

        node* n = 0;
        n = p->getFlow() > 1e-20 ? p->getNodeIn() : p->getNodeOut();
        if (n != 0) computeConvectiveStimuliRecursive(n, addedStimulus);
      }
    }
  }
}

void network::calculateConductedStimuli() {
  for (int j = 0; j < totalPores; ++j) {
    pore* p = getPore(j);
    if (!p->getClosed()) p->setConductedStim(0);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && abs(p->getFlow()) > 1e-20) {
      if (p->getConvectedStim() > 0) {
        for (int j = 0; j < totalPores; ++j) {
          pore* p = getPore(j);
          if (!p->getClosed()) p->setVisited(false);
        }

        double addedStimulus =
            log10(1 + (p->getConvectedStim() /
                       (abs(p->getFlow() * 6.0e13) + Qref * 6.0e13)));
        p->setConductedStim(p->getConductedStim() + addedStimulus);
        p->setVisited(true);

        node* n = 0;
        n = p->getFlow() > 1e-20 ? p->getNodeOut() : p->getNodeIn();
        if (n != 0) computeConductiveStimuliRecursive(n, addedStimulus);
      }
    }
  }
}

void network::computeConvectiveStimuliRecursive(node* n, double stimulus) {
  vector<int>& connectedPores = n->getConnectedPores();
  for (unsigned j = 0; j < connectedPores.size(); ++j) {
    pore* pp = getPore(connectedPores[j] - 1);
    if (!pp->getClosed() && !pp->getVisited()) {
      if ((pp->getNodeOut() == n && pp->getFlow() > 1e-20) ||
          (pp->getNodeIn() == n && pp->getFlow() < -1e-20)) {
        double newAddedStimulus =
            abs(pp->getFlow()) / n->getFlow() * stimulus * exp(-decayConv);

        pp->setConvectedStim(pp->getConvectedStim() + newAddedStimulus);
        pp->setVisited(true);

        node* nn = 0;
        nn = pp->getFlow() > 1e-20 ? pp->getNodeIn() : pp->getNodeOut();
        if (nn != 0 && newAddedStimulus > 0.01)
          computeConvectiveStimuliRecursive(nn, newAddedStimulus);
      }
    }
  }
}

void network::computeConductiveStimuliRecursive(node* n, double stimulus) {
  vector<int>& connectedPores = n->getConnectedPores();
  for (unsigned j = 0; j < connectedPores.size(); ++j) {
    pore* pp = getPore(connectedPores[j] - 1);
    if (!pp->getClosed() && !pp->getVisited()) {
      if ((pp->getNodeIn() == n && pp->getFlow() > 1e-20) ||
          (pp->getNodeOut() == n && pp->getFlow() < -1e-20)) {
        double newAddedStimulus = 1. / double(n->getFeedingVesselsNumber()) *
                                  stimulus * exp(-decayCond);
        pp->setConductedStim(pp->getConductedStim() + newAddedStimulus);
        pp->setVisited(true);

        node* nn = 0;
        nn = pp->getFlow() > 1e-20 ? pp->getNodeOut() : pp->getNodeIn();
        if (nn != 0 && newAddedStimulus > 0.01)
          computeConductiveStimuliRecursive(nn, newAddedStimulus);
      }
    }
  }
}

bool network::solvePressureInAngioModel() {
  assignBloodViscosities();
  assignVolumes();
  assignConductivities();

  // set constant flow rate
  double Q1(0), Q2(0), A, B;

  pressureIn = 1;
  pressureOut = 0;
  solvePressures();
  updateFlows();
  Q1 = getOutletFlow();

  pressureIn = 2;
  pressureOut = 0;
  solvePressures();
  updateFlows();
  Q2 = getOutletFlow();

  B = (Q1 - Q2 * 1 / 2) / (1 - 1 / 2);
  A = (Q1 - B);

  bool remodel = false;
  double oldDelta = deltaP;
  deltaP = flowRate / A - B / A;

  if (abs(deltaP - oldDelta) > 0.1) remodel = true;

  pressureIn = deltaP;
  pressureOut = 0;
  solvePressures();
  updateFlows();

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      vector<int> connectedPores = n->getConnectedPores();
      double sumFlowIn = 0;
      double massIn = 0;
      int vesselIn = 0;
      for (unsigned j = 0; j < connectedPores.size(); ++j) {
        pore* p = getPore(connectedPores[j] - 1);
        if (!p->getClosed()) {
          if ((p->getNodeIn() == n && p->getFlow() > 0) ||
              (p->getNodeOut() == n && p->getFlow() < 0)) {
            sumFlowIn += abs(p->getFlow());
            massIn += p->getHDConcentration() * abs(p->getFlow());
            vesselIn++;
          }
        }
      }
      n->setFeedingVesselsNumber(vesselIn);
      n->setFlow(sumFlowIn);
      n->setMassFlow(massIn);
    }
  }

  return remodel;
}

bool network::solvePressureInAngioModelWthPhaseSeparation() {
  bool remodel = false;
  double oldDelta = deltaP;

  assignBloodViscosities();
  assignVolumes();
  assignConductivities();

  // Solve pressure distribution
  bool stillMorePoresToClose = true;

  while (stillMorePoresToClose) {
    clusterEverything();

    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed() && p->getExist() == 't') {
        if (p->getClusterExist()->getSpanning() == false) {
          p->setConductivity(1e-200);
          p->setExist('f');
        }
      }
    }

    stillMorePoresToClose = false;

    // set constant flow rate
    double Q1(0), Q2(0), A, B;

    pressureIn = 1;
    pressureOut = 0;
    solvePressures();
    updateFlows();
    Q1 = getOutletFlow();

    pressureIn = 2;
    pressureOut = 0;
    solvePressures();
    updateFlows();
    Q2 = getOutletFlow();

    B = (Q1 - Q2 * 1 / 2) / (1 - 1 / 2);
    A = (Q1 - B);

    deltaP = flowRate / A - B / A;

    pressureIn = deltaP;
    pressureOut = 0;
    solvePressures();
    updateFlows();

    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      if (!n->getClosed()) {
        vector<int> connectedPores = n->getConnectedPores();
        double sumFlowIn = 0;
        double massIn = 0;
        int vesselIn = 0;
        for (unsigned j = 0; j < connectedPores.size(); ++j) {
          pore* p = getPore(connectedPores[j] - 1);
          if (!p->getClosed()) {
            if ((p->getNodeIn() == n && p->getFlow() > 1e-20) ||
                (p->getNodeOut() == n && p->getFlow() < -1e-20)) {
              sumFlowIn += abs(p->getFlow());
              massIn += p->getHDConcentration() * abs(p->getFlow());
              vesselIn++;
            }
          }
        }
        n->setFeedingVesselsNumber(vesselIn);
        n->setFlow(sumFlowIn);
        n->setMassFlow(massIn);
      }
    }

    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);

      if (!p->getClosed() && !p->getInlet() && !p->getOutlet() &&
          abs(p->getFlow()) > 1e-20) {
        node* n = 0;
        double FQB = 0;
        double FQE = 0;
        double Din = 0;
        double Dout = 0;
        double X0 = 0;
        double HDin = 0;
        double A(0), B(0);
        if (p->getFlow() > 1e-20) n = p->getNodeOut();
        if (p->getFlow() < -1e-20) n = p->getNodeIn();
        if (n != 0) {
          double numberOfFeedingVessels = 0;
          double numberOfFeededVessels = 0;

          if (n->getFlow() > 1e-20) FQB = abs(p->getFlow()) / n->getFlow();
          for (unsigned j = 0; j < n->getConnectedPores().size(); ++j) {
            pore* pp = getPore(n->getConnectedPores()[j] - 1);
            if (pp != p &&
                ((pp->getFlow() > 1e-20 && pp->getNodeIn() == n) ||
                 (pp->getFlow() < -1e-20 &&
                  pp->getNodeOut() == n)))  // pp belongs to a feeding vessel
            {
              Din += 2 * pp->getRadius();
              HDin += pp->getHDConcentration() * abs(pp->getFlow());
              numberOfFeedingVessels++;
            }
            if (pp != p &&
                ((pp->getFlow() > 1e-20 && pp->getNodeOut() == n) ||
                 (pp->getFlow() < -1e-20 &&
                  pp->getNodeIn() == n)))  // pp belongs to a feeded vessel
            {
              Dout += 2 * pp->getRadius();
              numberOfFeededVessels++;
            }
          }

          if (numberOfFeedingVessels != 0 && n->getFlow() > 1e-20) {
            Din /= double(numberOfFeedingVessels);
            HDin /= double(numberOfFeedingVessels) * n->getFlow();
            B = 1 + 6.98e-6 * (1 - HDin * 0.45) / Din;
          }

          if (numberOfFeededVessels != 0 && n->getFlow() > 1e-20)
            Dout /= double(numberOfFeededVessels);

          if (Din != 0 && Dout != 0) {
            X0 = 0.4e-6 / Din;
            A = -6.96e-6 * log(2 * p->getRadius() / Dout) / Din;
          }

          if (FQB <= X0)  // close pore momentarlily
          {
            p->setConductivity(1e-200);
            p->setExist('f');
            stillMorePoresToClose = true;
          } else {
            if (FQB + X0 > 0.9999) FQB = 0.9999 - X0;
            if (FQB > X0)
              FQE =
                  1 / (1 + exp(-(A + B * log((FQB - X0) / (1 - (FQB + X0))))));
          }
        }
        p->setFQE(FQE);
      } else
        p->setFQE(1);
    }

    // Thread management
    if (cancel) break;
  }

  if (abs(deltaP - oldDelta) > 0.1) remodel = true;

  // normalise FQE
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      vector<int> connectedPores = n->getConnectedPores();
      double sumFQE = 0;
      for (unsigned j = 0; j < connectedPores.size(); ++j) {
        pore* p = getPore(connectedPores[j] - 1);
        if (!p->getClosed()) {
          if ((p->getNodeOut() == n && p->getFlow() > 1e-20) ||
              (p->getNodeIn() == n && p->getFlow() < -1e-20))
            sumFQE += p->getFQE();
        }
      }
      if (sumFQE != 0)
        for (unsigned j = 0; j < connectedPores.size(); ++j) {
          pore* p = getPore(connectedPores[j] - 1);
          if (!p->getClosed()) {
            if ((p->getNodeOut() == n && p->getFlow() > 1e-20) ||
                (p->getNodeIn() == n && p->getFlow() < -1e-20)) {
              p->setFQE(p->getFQE() / sumFQE);
            }
          }
        }
    }
  }

  return remodel;
}
