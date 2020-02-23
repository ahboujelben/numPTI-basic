/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

using namespace std;

void network::runAngiogenesisOnLattice() {
  cout << "Starting Angiogenesis... " << endl;

  double startTime, endTime;
  startTime = tools::getCPUTime();

  // post-processing
  if (videoRecording) record = true;

  initialiseSimulation();
  setupTissueProperties();
  setupTAFDistribution();
  setupFNDistribution();
  assignInitialBloodViscosities();

  set<int> sproutTips;
  initialiseSroutTips(sproutTips);

  timeStep = 0.005;
  cout << "Time step (x 1.5day):  " << timeStep << endl;

  double timeSoFar = 0;
  double timeToRemodel = 0;

  while (timeSoFar < simulationTime) {
    timeSoFar += timeStep;
    timeToRemodel += timeStep;

    updateSproutTipPositions(sproutTips);  // update EC tips positions

    if (updateBlockAttributes)
      updateChemicalConcentrations();  // update TAF/FN/MDE concentrations in
                                       // tissue

    if (timeToRemodel > 1)  // flow and remodel vasculature radii every 1.5 day
    {
      remodelVasculature();
      timeToRemodel = 0;
    }

    setBranching(sproutTips);  // perform branching

    if (branchingWSS)
      setBranchingWSS(
          sproutTips);  // perform branching associated to shear stress

    emitPlotSignal();

    // Thread Management
    if (cancel) break;
  }

  // post-processing
  if (videoRecording) {
    record = false;
    extractVideo();
  }

  cout << "Simulation Time " << timeSoFar << endl;

  endTime = tools::getCPUTime();
  cout << "Processing Time: " << endTime - startTime << " s" << endl;
}

void network::setupTAFDistribution() {
  for (int i = 0; i < totalBlocks; ++i) {
    block *b = getBlock(i);
    if (!b->getClosed()) {
      double TAF(0);
      if (networkSource == 4)  // Tumour
      {
        if (circularTumour) {
          double mu = (sqrt(5) - 0.1) / (sqrt(5) - 1);
          double r =
              Nz < 5 ? sqrt(pow(b->getXCoordinate() / xEdgeLength - 1, 2) +
                            pow(b->getYCoordinate() / yEdgeLength - 0.5, 2))
                     : sqrt(pow(b->getXCoordinate() / xEdgeLength - 1, 2) +
                            pow(b->getYCoordinate() / yEdgeLength - 0.5, 2) +
                            pow(b->getZCoordinate() / zEdgeLength - 0.5, 2));
          if (r > 0 && r <= 0.1) TAF = 1;
          if (r > 0.1) TAF = pow((mu - r) / (mu - 0.1), 2);
        }
        if (linearTumour) {
          TAF = exp(-pow(1 - b->getXCoordinate() / xEdgeLength, 2) / 0.45);
        }
      }
      if (networkSource == 5)  // Retina
      {
        double mu = 1;
        double r = sqrt(pow(b->getXCoordinate() / xEdgeLength - 0.5, 2) +
                        pow(b->getYCoordinate() / yEdgeLength - 0.5, 2));
        TAF =
            1 - 0.45 * exp(-pow(2 * r, 2) / 0.45);  // exp(-pow(1-2*r,2)/0.45);
      }
      b->setTAFConcentration(TAF);
    }
  }
}

void network::setupFNDistribution() {
  for (int i = 0; i < totalBlocks; ++i) {
    block *b = getBlock(i);
    if (!b->getClosed()) {
      double FN(0);
      if (networkSource == 4)  // Tumour
      {
        FN = 0.75 * exp(-pow(b->getXCoordinate() / xEdgeLength, 2) / 0.45);
      }
      if (networkSource == 5)  // Retina
      {
        double mu = 1;
        double r = sqrt(pow(b->getXCoordinate() / xEdgeLength - 0.5, 2) +
                        pow(b->getYCoordinate() / yEdgeLength - 0.5, 2));
        FN = 0.1;  // 0.75*exp(-pow(2*r,2)/0.45);
      }
      b->setFNConcentration(FN);
    }
  }
}

void network::initialiseSroutTips(std::set<int> &sproutTips) {
  if (networkSource == 4)  // Tumour)
  {
    if (Nz < 5)
      for (int i = 1; i < initialTipsNumber + 1; ++i) {
        node *n = getNode(0, i * Ny / (initialTipsNumber + 1), 0);
        n->setSprouTip(true);
        sproutTips.insert(n->getId());
      }
    else
      for (int i = 1; i < initialTipsNumber + 1; ++i) {
        node *n = getNode(0, i * Ny / (initialTipsNumber + 1), Nz / 3);
        n->setSprouTip(true);
        sproutTips.insert(n->getId());
        n = getNode(0, i * Ny / (initialTipsNumber + 1), 2 * Nz / 3);
        n->setSprouTip(true);
        sproutTips.insert(n->getId());
      }
  }

  if (networkSource == 5)  // Retina
  {
    for (node *n : tableOfAllNodes) {
      if (n->getId() <= 2)  // gangliom
        continue;
      n->setSprouTip(true);
      sproutTips.insert(n->getId());
    }
  }
}

void network::calculateTimeStepForAngio() {
  // timeStep=min(min(min(1/(angio_Mu+4*angio_Epsilon/h_sq),1/angio_Eta),1/(angio_Gamma-angio_Beta)),h_sq/(4*(angio_D+Chi_func(0)+angio_Rho)));
  timeStep = 0.005;
  cout << "Time step (x 1.5day):  " << timeStep << endl;
}

void network::updateChemicalConcentrations() {
  double subTimeStep =
      min(min(min(1 / (angio_Mu + 4 * angio_Epsilon / h_sq), 1 / angio_Eta),
              1 / (angio_Gamma)),
          1 / angio_Beta) /
      2;

  double timeSoFar(0);
  while (timeSoFar < timeStep) {
    timeSoFar += subTimeStep;
    double dimensionPreFactor = Nz < 5 ? 4 : 6;
    // initialise maps

    vector<double> blockTAFConcentration;
    blockTAFConcentration.reserve(totalBlocks);
    vector<double> blockFNConcentration;
    blockFNConcentration.reserve(totalBlocks);
    vector<double> blockMDEConcentration;
    blockMDEConcentration.reserve(totalBlocks);
    double coeff1 = subTimeStep * angio_Alpha;
    double coeff2 = (1 - subTimeStep * angio_Mu -
                     dimensionPreFactor * subTimeStep * angio_Epsilon / h_sq);
    double coeff3 = subTimeStep * angio_Epsilon / h_sq;
    // update TAF, FN, MDE (maps)
    for (int i = 0; i < totalBlocks; ++i) {
      block *n = getBlock(i);
      if (!n->getClosed()) {
        double ii = n->getX();
        double jj = n->getY();
        double kk = n->getZ();

        double endothelialDensity = 0;

        node *s = getNode(ii, jj, kk);
        if (s != 0 && s->getSprouTip()) endothelialDensity = 1;

        double MDEconcentration = 0;

        if (endothelialDensity == 1) {
          blockTAFConcentration.push_back(n->getTAFConcentration() *
                                          (1 - subTimeStep * angio_Eta));
          blockFNConcentration.push_back(
              n->getFNConcentration() *
                  (1 - subTimeStep * angio_Gamma * n->getMDEConcentration()) +
              subTimeStep * angio_Beta);
          MDEconcentration = coeff1;
        }

        MDEconcentration += n->getMDEConcentration() * coeff2;

        block *nW, *nE, *nN, *nS, *nU, *nD;
        nW = getBlock(ii - 1, jj, kk);
        nE = getBlock(ii + 1, jj, kk);
        nS = getBlock(ii, jj - 1, kk);
        nN = getBlock(ii, jj + 1, kk);
        nD = getBlock(ii, jj, kk - 1);
        nU = getBlock(ii, jj, kk + 1);

        if (nW != 0) MDEconcentration += coeff3 * nW->getMDEConcentration();
        if (nE != 0) MDEconcentration += coeff3 * nE->getMDEConcentration();
        if (nS != 0) MDEconcentration += coeff3 * nS->getMDEConcentration();
        if (nN != 0) MDEconcentration += coeff3 * nN->getMDEConcentration();
        if (Nz >= 5)
          if (nD != 0) MDEconcentration += coeff3 * nD->getMDEConcentration();
        if (Nz >= 5)
          if (nU != 0) MDEconcentration += coeff3 * nU->getMDEConcentration();

        if (nW == 0) MDEconcentration += coeff3 * n->getMDEConcentration();
        if (nE == 0) MDEconcentration += coeff3 * n->getMDEConcentration();
        if (nS == 0) MDEconcentration += coeff3 * n->getMDEConcentration();
        if (nN == 0) MDEconcentration += coeff3 * n->getMDEConcentration();
        if (Nz >= 5)
          if (nD == 0) MDEconcentration += coeff3 * n->getMDEConcentration();
        if (Nz >= 5)
          if (nU == 0) MDEconcentration += coeff3 * n->getMDEConcentration();

        blockMDEConcentration.push_back(MDEconcentration);
      }
    }

    unsigned j(0);
    for (int i = 0; i < totalBlocks; ++i) {
      block *n = getBlock(i);
      if (!n->getClosed()) {
        n->setMDEConcentration(blockMDEConcentration[j]);
        j++;
        if (n->getMDEConcentration() < 0 || n->getMDEConcentration() > 1) {
          cout << "MDE concentration out of range: " << n->getMDEConcentration()
               << endl;
          cancel = true;
        }
      }
    }
  }
}

void network::updateSproutTipPositions(std::set<int> &sproutTips) {
  set<int> newSproutTips;

  for (set<int>::iterator it = sproutTips.begin(); it != sproutTips.end();
       ++it) {
    node *n = getNode(*it - 1);

    if (networkSource == 5)  // Retina
    {
      // prevent growth outside circular boundary
      double r_sq = pow(n->getXCoordinate() / xEdgeLength - 0.5, 2) +
                    pow(n->getYCoordinate() / yEdgeLength - 0.5, 2);
      if (r_sq > pow(0.5, 2)) continue;
    }

    // calculate P1..6
    vector<double> probabilities = generateEndothelialCellProbabilities(n);

    double sumProbabilities(0);
    for (int i = 0; i < 7; ++i) sumProbabilities += probabilities[i];

    double dice = uniform_real();
    double cumulutativeProbability = 0;
    int movementDirection = 0;
    for (int i = 0; i < 7; ++i) {
      cumulutativeProbability += probabilities[i] / sumProbabilities;
      if (dice < cumulutativeProbability) {
        movementDirection = i;
        break;
      }
    }

    // determine node to move to
    if (movementDirection == 0) {
      if (probabilities[7] > 1e-5) newSproutTips.insert(n->getId());
    }

    if (movementDirection == 2) {
      if (n->getIndexX() < Nx - 1 &&
          getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ()) == 0) {
        node *newN =
            addNode(n, n->getIndexX() + 1, n->getIndexY(), n->getIndexZ());
        newSproutTips.insert(newN->getId());

        addVessel(newN, n, 1);
      } else if (n->getIndexX() < Nx - 1 &&
                 getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ()) !=
                     0) {
        node *nn = getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ());
        int i = nn->getIndexX();
        int j = nn->getIndexY();
        int k = nn->getIndexZ();
        if (getPoreX(i, j, k) == 0) addVessel(nn, n, 1);
        newSproutTips.insert(nn->getId());
      }
    }

    if (movementDirection == 1) {
      if (n->getIndexX() > 0 &&
          getNode(n->getIndexX() - 1, n->getIndexY(), n->getIndexZ()) == 0) {
        node *newN =
            addNode(n, n->getIndexX() - 1, n->getIndexY(), n->getIndexZ());
        newSproutTips.insert(newN->getId());

        addVessel(n, newN, 1);
      } else if (n->getIndexX() > 0 &&
                 getNode(n->getIndexX() - 1, n->getIndexY(), n->getIndexZ()) !=
                     0) {
        node *nn = getNode(n->getIndexX() - 1, n->getIndexY(), n->getIndexZ());
        int i = n->getIndexX();
        int j = n->getIndexY();
        int k = n->getIndexZ();
        if (getPoreX(i, j, k) == 0) addVessel(n, nn, 1);
        newSproutTips.insert(nn->getId());
      }
    }

    if (movementDirection == 4) {
      if (n->getIndexY() < Ny - 1 &&
          getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ()) == 0) {
        node *newN =
            addNode(n, n->getIndexX(), n->getIndexY() + 1, n->getIndexZ());
        newSproutTips.insert(newN->getId());

        addVessel(newN, n, 2);
      } else if (n->getIndexY() < Ny - 1 &&
                 getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ()) !=
                     0) {
        node *nn = getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ());
        int i = nn->getIndexX();
        int j = nn->getIndexY();
        int k = nn->getIndexZ();
        if (getPoreY(i, j, k) == 0) addVessel(nn, n, 2);
        newSproutTips.insert(nn->getId());
      }
    }

    if (movementDirection == 3) {
      if (n->getIndexY() > 0 &&
          getNode(n->getIndexX(), n->getIndexY() - 1, n->getIndexZ()) == 0) {
        node *newN =
            addNode(n, n->getIndexX(), n->getIndexY() - 1, n->getIndexZ());
        newSproutTips.insert(newN->getId());

        addVessel(n, newN, 2);
      } else if (n->getIndexY() > 0 &&
                 getNode(n->getIndexX(), n->getIndexY() - 1, n->getIndexZ()) !=
                     0) {
        node *nn = getNode(n->getIndexX(), n->getIndexY() - 1, n->getIndexZ());
        int i = n->getIndexX();
        int j = n->getIndexY();
        int k = n->getIndexZ();
        if (getPoreY(i, j, k) == 0) addVessel(n, nn, 2);
        newSproutTips.insert(nn->getId());
      }
    }

    if (movementDirection == 6) {
      if (n->getIndexZ() < Nz - 1 &&
          getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1) == 0) {
        node *newN =
            addNode(n, n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1);
        newSproutTips.insert(newN->getId());

        addVessel(newN, n, 3);
      } else if (n->getIndexZ() < Nz - 1 &&
                 getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1) !=
                     0) {
        node *nn = getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1);
        int i = nn->getIndexX();
        int j = nn->getIndexY();
        int k = nn->getIndexZ();
        if (getPoreZ(i, j, k) == 0) addVessel(nn, n, 3);
        newSproutTips.insert(nn->getId());
      }
    }

    if (movementDirection == 5) {
      if (n->getIndexZ() > 0 &&
          getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1) == 0) {
        node *newN =
            addNode(n, n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1);
        newSproutTips.insert(newN->getId());

        addVessel(n, newN, 3);
      } else if (n->getIndexZ() > 0 && getNode(n->getIndexX(), n->getIndexY(),
                                               n->getIndexZ() - 1) != 0) {
        node *nn = getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1);
        int i = n->getIndexX();
        int j = n->getIndexY();
        int k = n->getIndexZ();
        if (getPoreZ(i, j, k) == 0) addVessel(n, nn, 3);
        newSproutTips.insert(nn->getId());
      }
    }
  }

  for (set<int>::iterator it = sproutTips.begin(); it != sproutTips.end();
       ++it) {
    node *n = getNode(*it - 1);
    n->setSprouTip(false);
  }

  for (set<int>::iterator it = newSproutTips.begin(); it != newSproutTips.end();
       ++it) {
    node *n = getNode(*it - 1);
    n->setSprouTip(true);
    n->setAge(n->getAge() + timeStep);
  }

  sproutTips = newSproutTips;

  updateRanking();
}

void network::setBranching(std::set<int> &sproutTips) {
  for (set<int>::iterator it = sproutTips.begin(); it != sproutTips.end();
       ++it) {
    node *n = getNode(*it - 1);
    block *b = getBlock(n->getIndexX(), n->getIndexY(), n->getIndexZ());

    if (n->getAge() > angio_Psi) {
      double branchingProbability(0);
      if (b->getTAFConcentration() > 0.3 && b->getTAFConcentration() <= 0.5)
        branchingProbability = 0.2;
      if (b->getTAFConcentration() > 0.5 && b->getTAFConcentration() <= 0.7)
        branchingProbability = 0.3;
      if (b->getTAFConcentration() > 0.7 && b->getTAFConcentration() <= 0.8)
        branchingProbability = 0.4;
      if (b->getTAFConcentration() > 0.8) branchingProbability = 1;

      double dice = uniform_real();
      if (dice < branchingProbability) {
        // starting braching
        generateEndothelialCellProbabilities(n);
        vector<double> probabilities = generateEndothelialCellProbabilities(n);

        double sumProbabilities(0);
        for (int i = 1; i < 7; ++i) sumProbabilities += probabilities[i];

        double dice = uniform_real();
        double cumulutativeProbability = 0;
        int movementDirection = 0;
        for (int i = 1; i < 5; ++i) {
          cumulutativeProbability += probabilities[i] / sumProbabilities;
          if (dice < cumulutativeProbability) {
            movementDirection = i;
            break;
          }
        }

        if (movementDirection == 2) {
          if (n->getIndexX() < Nx - 1 &&
              getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ()) ==
                  0) {
            node *newN =
                addNode(n, n->getIndexX() + 1, n->getIndexY(), n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(newN, n, 1);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexX() < Nx - 1 &&
                     getNode(n->getIndexX() + 1, n->getIndexY(),
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ());
            // sproutTips.insert(nn);

            int i = nn->getIndexX();
            int j = nn->getIndexY();
            int k = nn->getIndexZ();
            if (getPoreX(i, j, k) == 0) addVessel(nn, n, 1);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 1) {
          if (n->getIndexX() > 0 && getNode(n->getIndexX() - 1, n->getIndexY(),
                                            n->getIndexZ()) == 0) {
            node *newN =
                addNode(n, n->getIndexX() - 1, n->getIndexY(), n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(n, newN, 1);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexX() > 0 &&
                     getNode(n->getIndexX() - 1, n->getIndexY(),
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX() - 1, n->getIndexY(), n->getIndexZ());
            // sproutTips.insert(nn);

            int i = n->getIndexX();
            int j = n->getIndexY();
            int k = n->getIndexZ();
            if (getPoreX(i, j, k) == 0) addVessel(n, nn, 1);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 4) {
          if (n->getIndexY() < Ny - 1 &&
              getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ()) ==
                  0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY() + 1, n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(newN, n, 2);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexY() < Ny - 1 &&
                     getNode(n->getIndexX(), n->getIndexY() + 1,
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ());
            // sproutTips.insert(nn);

            int i = nn->getIndexX();
            int j = nn->getIndexY();
            int k = nn->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(nn, n, 2);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 3) {
          if (n->getIndexY() > 0 && getNode(n->getIndexX(), n->getIndexY() - 1,
                                            n->getIndexZ()) == 0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY() - 1, n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(n, newN, 2);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexY() > 0 &&
                     getNode(n->getIndexX(), n->getIndexY() - 1,
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY() - 1, n->getIndexZ());
            // sproutTips.insert(nn);

            int i = n->getIndexX();
            int j = n->getIndexY();
            int k = n->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(n, nn, 2);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 6) {
          if (n->getIndexZ() < Nz - 1 && getNode(n->getIndexX(), n->getIndexY(),
                                                 n->getIndexZ() + 1) == 0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1);
            sproutTips.insert(newN->getId());

            addVessel(newN, n, 3);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexZ() < Nz - 1 &&
                     getNode(n->getIndexX(), n->getIndexY(),
                             n->getIndexZ() + 1) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1);
            // sproutTips.insert(nn);

            int i = nn->getIndexX();
            int j = nn->getIndexY();
            int k = nn->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(nn, n, 3);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 5) {
          if (n->getIndexZ() > 0 && getNode(n->getIndexX(), n->getIndexY(),
                                            n->getIndexZ() - 1) == 0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1);
            sproutTips.insert(newN->getId());

            addVessel(n, newN, 3);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexZ() > 0 &&
                     getNode(n->getIndexX(), n->getIndexY(),
                             n->getIndexZ() - 1) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1);
            // sproutTips.insert(nn);

            int i = n->getIndexX();
            int j = n->getIndexY();
            int k = n->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(n, nn, 3);

            n->setAge(0);
            nn->setAge(0);
          }
        }
      }
    }
  }
  updateRanking();
}

void network::setBranchingWSS(std::set<int> &sproutTips) {
  for (int i = 0; i < totalPores; ++i) {
    pore *p = getPore(i);

    if (!p->getClosed() && !p->getParentVessel() && abs(p->getFlow()) > 1e-20) {
      node *n = 0;
      if (p->getFlow() > 0)
        n = p->getNodeIn();
      else
        n = p->getNodeOut();
      block *b = getBlock(n->getIndexX(), n->getIndexY(), n->getIndexZ());

      double WSSRatio = p->getWSS() / angio_TauMax;
      double branchingProbability(0);
      if (b->getTAFConcentration() > 0.3 && b->getTAFConcentration() <= 0.5) {
        if (WSSRatio > 0.2 && WSSRatio <= 0.4) branchingProbability = 0.02;
        if (WSSRatio > 0.4 && WSSRatio <= 0.6) branchingProbability = 0.04;
        if (WSSRatio > 0.6 && WSSRatio <= 0.8) branchingProbability = 0.06;
        if (WSSRatio > 0.8 && WSSRatio <= 1.0) branchingProbability = 0.08;
      }
      if (b->getTAFConcentration() > 0.5 && b->getTAFConcentration() <= 0.7) {
        if (WSSRatio > 0.2 && WSSRatio <= 0.4) branchingProbability = 0.03;
        if (WSSRatio > 0.4 && WSSRatio <= 0.6) branchingProbability = 0.06;
        if (WSSRatio > 0.6 && WSSRatio <= 0.8) branchingProbability = 0.09;
        if (WSSRatio > 0.8 && WSSRatio <= 1.0) branchingProbability = 0.12;
      }
      if (b->getTAFConcentration() > 0.7 && b->getTAFConcentration() <= 0.8) {
        if (WSSRatio > 0.2 && WSSRatio <= 0.4) branchingProbability = 0.04;
        if (WSSRatio > 0.4 && WSSRatio <= 0.6) branchingProbability = 0.08;
        if (WSSRatio > 0.6 && WSSRatio <= 0.8) branchingProbability = 0.12;
        if (WSSRatio > 0.8 && WSSRatio <= 1.0) branchingProbability = 0.16;
      }
      if (b->getTAFConcentration() > 0.8) {
        if (WSSRatio > 0.2 && WSSRatio <= 0.4) branchingProbability = 0.10;
        if (WSSRatio > 0.4 && WSSRatio <= 0.6) branchingProbability = 0.20;
        if (WSSRatio > 0.6 && WSSRatio <= 0.8) branchingProbability = 0.30;
        if (WSSRatio > 0.8 && WSSRatio <= 1.0) branchingProbability = 0.40;
      }

      double dice = uniform_real();
      if (dice < branchingProbability) {
        // starting braching
        generateEndothelialCellProbabilities(n);
        vector<double> probabilities = generateEndothelialCellProbabilities(n);

        double sumProbabilities(0);
        for (int i = 1; i < 7; ++i) sumProbabilities += probabilities[i];

        double dice = uniform_real();
        double cumulutativeProbability = 0;
        int movementDirection = 0;
        for (int i = 1; i < 5; ++i) {
          cumulutativeProbability += probabilities[i] / sumProbabilities;
          if (dice < cumulutativeProbability) {
            movementDirection = i;
            break;
          }
        }

        if (movementDirection == 2) {
          if (n->getIndexX() < Nx - 1 &&
              getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ()) ==
                  0) {
            node *newN =
                addNode(n, n->getIndexX() + 1, n->getIndexY(), n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(newN, n, 1);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexX() < Nx - 1 &&
                     getNode(n->getIndexX() + 1, n->getIndexY(),
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX() + 1, n->getIndexY(), n->getIndexZ());
            // sproutTips.insert(nn);

            int i = nn->getIndexX();
            int j = nn->getIndexY();
            int k = nn->getIndexZ();
            if (getPoreX(i, j, k) == 0) addVessel(nn, n, 1);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 1) {
          if (n->getIndexX() > 0 && getNode(n->getIndexX() - 1, n->getIndexY(),
                                            n->getIndexZ()) == 0) {
            node *newN =
                addNode(n, n->getIndexX() - 1, n->getIndexY(), n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(n, newN, 1);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexX() > 0 &&
                     getNode(n->getIndexX() - 1, n->getIndexY(),
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX() - 1, n->getIndexY(), n->getIndexZ());
            // sproutTips.insert(nn);

            int i = n->getIndexX();
            int j = n->getIndexY();
            int k = n->getIndexZ();
            if (getPoreX(i, j, k) == 0) addVessel(n, nn, 1);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 4) {
          if (n->getIndexY() < Ny - 1 &&
              getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ()) ==
                  0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY() + 1, n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(newN, n, 2);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexY() < Ny - 1 &&
                     getNode(n->getIndexX(), n->getIndexY() + 1,
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY() + 1, n->getIndexZ());
            // sproutTips.insert(nn);

            int i = nn->getIndexX();
            int j = nn->getIndexY();
            int k = nn->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(nn, n, 2);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 3) {
          if (n->getIndexY() > 0 && getNode(n->getIndexX(), n->getIndexY() - 1,
                                            n->getIndexZ()) == 0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY() - 1, n->getIndexZ());
            sproutTips.insert(newN->getId());

            addVessel(n, newN, 2);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexY() > 0 &&
                     getNode(n->getIndexX(), n->getIndexY() - 1,
                             n->getIndexZ()) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY() - 1, n->getIndexZ());
            // sproutTips.insert(nn);

            int i = n->getIndexX();
            int j = n->getIndexY();
            int k = n->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(n, nn, 2);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 6) {
          if (n->getIndexZ() < Nz - 1 && getNode(n->getIndexX(), n->getIndexY(),
                                                 n->getIndexZ() + 1) == 0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1);
            sproutTips.insert(newN->getId());

            addVessel(newN, n, 3);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexZ() < Nz - 1 &&
                     getNode(n->getIndexX(), n->getIndexY(),
                             n->getIndexZ() + 1) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() + 1);
            // sproutTips.insert(nn);

            int i = nn->getIndexX();
            int j = nn->getIndexY();
            int k = nn->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(nn, n, 3);

            n->setAge(0);
            nn->setAge(0);
          }
        }

        if (movementDirection == 5) {
          if (n->getIndexZ() > 0 && getNode(n->getIndexX(), n->getIndexY(),
                                            n->getIndexZ() - 1) == 0) {
            node *newN =
                addNode(n, n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1);
            sproutTips.insert(newN->getId());

            addVessel(n, newN, 3);

            n->setAge(0);
            newN->setAge(0);
          } else if (n->getIndexZ() > 0 &&
                     getNode(n->getIndexX(), n->getIndexY(),
                             n->getIndexZ() - 1) != 0) {
            node *nn =
                getNode(n->getIndexX(), n->getIndexY(), n->getIndexZ() - 1);
            // sproutTips.insert(nn);

            int i = n->getIndexX();
            int j = n->getIndexY();
            int k = n->getIndexZ();
            if (getPoreY(i, j, k) == 0) addVessel(n, nn, 3);

            n->setAge(0);
            nn->setAge(0);
          }
        }
      }
    }
  }
  updateRanking();
}

void network::remodelVasculature() {
  bool unstablePressureField = true;
  double PVsInjectedSoFar(0);
  double dampingFactor(1);  // used to help radii adaptation converge

  // solve pressure in vasculature
  phaseSeparation ? solvePressureInAngioModelWthPhaseSeparation()
                  : solvePressureInAngioModel();

  while (unstablePressureField) {
    double volumeInjectedSoFar = 0;
    while (volumeInjectedSoFar / totalPoresVolume <
           1)  // inject 2 PV before recalculating radii
    {
      volumeInjectedSoFar += phaseSeparation
                                 ? runHaematocritFlowWithPhaseSeparation()
                                 : runHaematocritFlow();
      // Thread management
      if (cancel) break;
    }

    PVsInjectedSoFar += volumeInjectedSoFar / totalPoresVolume;

    if (shuntPrevention) {
      calculateConvectedStimuli();
      calculateConductedStimuli();
    }

    bool unstableVascularStructure = true;
    while (unstableVascularStructure) {
      // adjust vascular radii until a stable structure is reached
      unstableVascularStructure = recalculateRadii(dampingFactor);

      // Thread management
      if (cancel) break;
    }
    dampingFactor += 0.1;

    unstablePressureField = phaseSeparation
                                ? solvePressureInAngioModelWthPhaseSeparation()
                                : solvePressureInAngioModel();

    // A safety control to prevent non-convergence
    if (PVsInjectedSoFar > 50) break;

    // Thread Management
    if (cancel) break;
  }
}

bool network::recalculateRadii(double it) {
  bool recalculate = false;
  for (int i = 0; i < totalPores; ++i) {
    pore *p = getPore(i);
    if (!p->getClosed() && !p->getParentVessel() && abs(p->getFlow()) > 1e-20 &&
        p->getHDConcentration() > 0.05) {
      // WSS
      double tauW = 4 * p->getViscosity() /
                    (tools::pi() * pow(p->getRadius(), 3)) * abs(p->getFlow());
      double Swss = log10(10 * tauW + tauRef);
      p->setWSS(tauW);
      // cout<<p->getFlow()<<" "<<wss<<" "<<log10(10*wss+0.103)<<endl;

      // Intravascular pressure
      double tauE =
          0.1 *
          (100.0 -
           86.0 * exp(-5000.0 * pow(log10(log10(max(
                                        p->getAveragePressure() / 133, 10.1))),
                                    5.4)));
      double Sp = -Kp * log10(10 * tauE);
      // cout<<tauE<<" "<<p->getAveragePressure()<<"
      // "<<-0.1*log10(10*tauE)<<endl;

      // Metabolic stimulus
      double Sm = 0;

      if (!shuntPrevention)
        Sm = p->getHDConcentration() > 0.001
                 ? Km * log10((flowRate / (abs(p->getFlow()) *
                                           p->getHDConcentration() * 0.45)) +
                              1.0)
                 : 0;
      else {
        double SmConv = p->getConvectedStim() > 0
                            ? Km * log10(1 + p->getConvectedStim() /
                                                 (abs(p->getFlow() * 6.0e13) +
                                                  Qref * 6.0e13))
                            : 0;
        // cout<<"Smconv "<<SmConv<<endl;
        double SmCond =
            Km * Kc * (p->getConductedStim() / (p->getConductedStim() + J0));
        // cout<<"Smcond "<<SmCond<<endl;
        // cout<<SmConv<<" "<<SmCond<<endl;
        Sm = SmCond + SmConv;
      }

      double deltaR =
          it == 0 ? (Swss + Sp + Sm - Ks) * p->getRadius() * timeStep
                  : (Swss + Sp + Sm - Ks) * p->getRadius() * timeStep / it;
      // cout<<Swss<<" "<<Sp<<" "<<Sm<<" "<<endl;
      double newRadius = p->getRadius() + deltaR;
      if (newRadius > 12e-6) newRadius = 12e-6;
      if (newRadius < 2e-6) newRadius = 2e-6;

      if (abs(p->getRadius() - newRadius) > 0.001e-6) recalculate = true;

      p->setRadius(newRadius);
    }
  }
  return recalculate;
}

double network::Chi_func(double c) { return angio_Chi / (1 + angio_Delta * c); }

std::vector<double> network::generateEndothelialCellProbabilities(node *no) {
  double ii = no->getIndexX();
  double jj = no->getIndexY();
  double kk = no->getIndexZ();

  block *n = getBlock((ii, jj, kk));
  block *nW, *nE, *nN, *nS, *nU, *nD;
  nW = getBlock(ii - 1, jj, kk);
  nE = getBlock(ii + 1, jj, kk);
  nS = getBlock(ii, jj - 1, kk);
  nN = getBlock(ii, jj + 1, kk);
  nD = getBlock(ii, jj, kk - 1);
  nU = getBlock(ii, jj, kk + 1);

  double dimensionPreFactor = Nz < 5 ? 4 : 6;
  vector<double> probabilities;

  double c = n->getTAFConcentration();
  double f = n->getFNConcentration();
  double cW = nW == 0 ? n->getTAFConcentration() : nW->getTAFConcentration();
  double fW = nW == 0 ? n->getFNConcentration() : nW->getFNConcentration();
  double cE = nE == 0 ? n->getTAFConcentration() : nE->getTAFConcentration();
  double fE = nE == 0 ? n->getFNConcentration() : nE->getFNConcentration();
  double cS = nS == 0 ? n->getTAFConcentration() : nS->getTAFConcentration();
  double fS = nS == 0 ? n->getFNConcentration() : nS->getFNConcentration();
  double cN = nN == 0 ? n->getTAFConcentration() : nN->getTAFConcentration();
  double fN = nN == 0 ? n->getFNConcentration() : nN->getFNConcentration();
  double cD = nD == 0 ? n->getTAFConcentration() : nD->getTAFConcentration();
  double fD = nD == 0 ? n->getFNConcentration() : nD->getFNConcentration();
  double cU = nU == 0 ? n->getTAFConcentration() : nU->getTAFConcentration();
  double fU = nU == 0 ? n->getFNConcentration() : nU->getFNConcentration();

  double P0 =
      1 -
      dimensionPreFactor * timeStep * angio_D /
          h_sq;  //+timeStep*angio_Alpha*Chi_func(c)/(4*h_sq*(1+angio_Alpha*c))*(pow(cE-cW,2)+pow(cN-cS,2))
                 //-timeStep*Chi_func(c)/h_sq*(cE+cW+cN+cS-4*c)-timeStep*angio_Rho/h_sq*(fE+fW+fN+fS-4*f);
  double P1 = timeStep * angio_D / h_sq -
              timeStep / (dimensionPreFactor * h_sq) *
                  (Chi_func(c) * (cE - cW) + angio_Rho * (fE - fW));
  double P2 = timeStep * angio_D / h_sq +
              timeStep / (dimensionPreFactor * h_sq) *
                  (Chi_func(c) * (cE - cW) + angio_Rho * (fE - fW));
  double P3 = timeStep * angio_D / h_sq -
              timeStep / (dimensionPreFactor * h_sq) *
                  (Chi_func(c) * (cN - cS) + angio_Rho * (fN - fS));
  double P4 = timeStep * angio_D / h_sq +
              timeStep / (dimensionPreFactor * h_sq) *
                  (Chi_func(c) * (cN - cS) + angio_Rho * (fN - fS));
  double P5 = Nz < 5
                  ? 0
                  : timeStep * angio_D / h_sq -
                        timeStep / (dimensionPreFactor * h_sq) *
                            (Chi_func(c) * (cU - cD) + angio_Rho * (fU - fD));
  double P6 = Nz < 5
                  ? 0
                  : timeStep * angio_D / h_sq +
                        timeStep / (dimensionPreFactor * h_sq) *
                            (Chi_func(c) * (cU - cD) + angio_Rho * (fU - fD));

  // avoid hitting its own sprout
  if (getNode(ii - 1, jj, kk) != 0 && getPoreX(ii, jj, kk) != 0 ||
      no->getIndexX() == 0)
    P1 = 0;
  if (getNode(ii + 1, jj, kk) != 0 && getPoreX(ii + 1, jj, kk) != 0 ||
      no->getIndexX() == Nx - 1)
    P2 = 0;
  if (getNode(ii, jj - 1, kk) != 0 && getPoreY(ii, jj, kk) != 0 ||
      no->getIndexX() == 0 || no->getIndexY() == 0)
    P3 = 0;
  if (getNode(ii, jj + 1, kk) != 0 && getPoreY(ii, jj + 1, kk) != 0 ||
      no->getIndexX() == 0 || no->getIndexY() == Ny - 1)
    P4 = 0;
  if (getNode(ii, jj, kk - 1) != 0 && getPoreZ(ii, jj, kk) != 0 ||
      no->getIndexX() == 0 || no->getIndexZ() == 0)
    P5 = 0;
  if (getNode(ii, jj, kk + 1) != 0 && getPoreZ(ii, jj, kk + 1) != 0 ||
      no->getIndexX() == 0 || no->getIndexZ() == Nz - 1)
    P6 = 0;

  double P7(0);

  // cout<<P1<<" "<<P2<<" "<<P3<<" "<<P4<<endl;

  double add1(0), add2(0), add3(0), add4(0), add5(0), add6(0);
  if (P1 < 0) {
    add1 = -P1;
  }
  if (P2 < 0) {
    add2 = -P2;
  }
  if (P3 < 0) {
    add3 = -P3;
  }
  if (P4 < 0) {
    add4 = -P4;
  }
  if (P5 < 0) {
    add5 = -P5;
  }
  if (P6 < 0) {
    add6 = -P6;
  }

  if (P1 > 0) P1 += add2;
  if (P2 > 0) P2 += add1;
  if (P3 > 0) P3 += add4;
  if (P4 > 0) P4 += add3;
  if (P5 > 0) P5 += add6;
  if (P6 > 0) P6 += add5;

  if (P0 < 0) {
    P0 = 0;
  }
  if (P1 < 0) {
    P1 = 0;
  }
  if (P2 < 0) {
    P2 = 0;
  }
  if (P3 < 0) {
    P3 = 0;
  }
  if (P4 < 0) {
    P4 = 0;
  }
  if (P5 < 0) {
    P5 = 0;
  }
  if (P6 < 0) {
    P6 = 0;
  }

  double cumulutativeProbability = P0 + P1 + P2 + P3 + P4 + P5 + P6;
  if (cumulutativeProbability != 0) {
    P0 = P0 / cumulutativeProbability;
    P1 = P1 / cumulutativeProbability;
    P2 = P2 / cumulutativeProbability;
    P3 = P3 / cumulutativeProbability;
    P4 = P4 / cumulutativeProbability;
    P5 = P5 / cumulutativeProbability;
    P6 = P6 / cumulutativeProbability;
  }
  P7 = P1 + P2 + P3 + P4 + P5 + P6;

  probabilities.push_back(P0);
  probabilities.push_back(P1);
  probabilities.push_back(P2);
  probabilities.push_back(P3);
  probabilities.push_back(P4);
  probabilities.push_back(P5);
  probabilities.push_back(P6);
  probabilities.push_back(P7);

  return probabilities;
}

node *network::addNode(node *source, int i, int j, int k) {
  node *n = new node(i, j, k);
  n->setId(totalNodes + 1);
  n->setXCoordinate(n->getIndexX() * length);
  n->setYCoordinate(n->getIndexY() * length);
  n->setZCoordinate(n->getIndexZ() * length);

  // distort
  n->setXCoordinate(n->getXCoordinate() +
                    length * degreeOfDistortion * (-1 + 2 * uniform_real()));
  n->setYCoordinate(n->getYCoordinate() +
                    length * degreeOfDistortion * (-1 + 2 * uniform_real()));
  if (Nz != 1)
    n->setZCoordinate(n->getZCoordinate() +
                      length * degreeOfDistortion * (-1 + 2 * uniform_real()));

  tableOfAllNodes.push_back(n);
  tableOfNodes[i][j][k] = n;
  totalNodes++;
  totalOpenedNodes++;

  // set neighbours
  n->getNeighboors().push_back(source->getId());
  source->getNeighboors().push_back(n->getId());

  // inhertit the same age of the parent sprout
  n->setAge(source->getAge());

  return n;
}

pore *network::addVessel(node *nodeIn, node *nodeOut, int i) {
  pore *p = new pore(nodeIn, nodeOut);
  tableOfAllPores.push_back(p);
  p->setId(totalPores + 1);
  totalPores++;
  totalOpenedPores++;

  if (i == 1)  // X pore
    tableOfPoresX[nodeIn->getIndexX()][nodeIn->getIndexY()]
                 [nodeIn->getIndexZ()] = p;

  if (i == 2)  // Y pore
    tableOfPoresY[nodeIn->getIndexX()][nodeIn->getIndexY()]
                 [nodeIn->getIndexZ()] = p;

  if (i == 3)  // Z pore
    tableOfPoresZ[nodeIn->getIndexX()][nodeIn->getIndexY()]
                 [nodeIn->getIndexZ()] = p;

  double averageRadius(0);
  // set neighboors
  nodeIn->getConnectedPores().push_back(p->getId());
  nodeOut->getConnectedPores().push_back(p->getId());

  if (find(nodeIn->getNeighboors().begin(), nodeIn->getNeighboors().end(),
           nodeOut->getId()) == nodeIn->getNeighboors().end())
    nodeIn->getNeighboors().push_back(nodeOut->getId());
  if (find(nodeOut->getNeighboors().begin(), nodeOut->getNeighboors().end(),
           nodeIn->getId()) == nodeOut->getNeighboors().end())
    nodeOut->getNeighboors().push_back(nodeIn->getId());

  for (int i = 0; i < nodeIn->getConnectedPores().size(); ++i) {
    pore *pp = getPore(nodeIn->getConnectedPores()[i] - 1);
    if (pp != p) {
      p->getNeighboors().push_back(pp);
      pp->getNeighboors().push_back(p);
      averageRadius += pp->getRadius();
    }
  }

  for (int i = 0; i < nodeOut->getConnectedPores().size(); ++i) {
    pore *pp = getPore(nodeOut->getConnectedPores()[i] - 1);
    if (pp != p) {
      p->getNeighboors().push_back(pp);
      pp->getNeighboors().push_back(p);
      averageRadius += pp->getRadius();
    }
  }

  averageRadius /= double(p->getNeighboors().size());

  p->setRadius(min(6e-6, averageRadius));

  double length = sqrt(
      pow(p->getNodeIn()->getXCoordinate() - p->getNodeOut()->getXCoordinate(),
          2) +
      pow(p->getNodeIn()->getYCoordinate() - p->getNodeOut()->getYCoordinate(),
          2) +
      pow(p->getNodeIn()->getZCoordinate() - p->getNodeOut()->getZCoordinate(),
          2));
  p->setFullLength(length);
  p->setLength(length);

  p->setMembranePermeability(PVT);

  // Update tissue connection
  if (nodeIn != 0) {
    block *b1 =
        getBlock(nodeIn->getIndexX(), nodeIn->getIndexY(), nodeIn->getIndexZ());
    b1->neighbooringVesselsArea()[p->getId()] =
        tools::pi() * p->getRadius() * p->getLength();
    b1->setVolume(b1->getVolume() - p->getVolume() / 2);
    b1->setConnectedToVessel(true);
    p->neighbooringBlocksArea()[b1->getId()] =
        tools::pi() * p->getRadius() * p->getLength();
  }
  if (nodeOut != 0) {
    block *b2 = getBlock(nodeOut->getIndexX(), nodeOut->getIndexY(),
                         nodeOut->getIndexZ());
    b2->neighbooringVesselsArea()[p->getId()] =
        tools::pi() * p->getRadius() * p->getLength();
    b2->setVolume(b2->getVolume() - p->getVolume() / 2);
    b2->setConnectedToVessel(true);
    p->neighbooringBlocksArea()[b2->getId()] =
        tools::pi() * p->getRadius() * p->getLength();
  }

  return p;
}
