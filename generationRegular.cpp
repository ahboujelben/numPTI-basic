/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

using namespace std;

void network::setupRegularModel() {
  cout << "Creating Nodes..." << endl;
  createNodes();
  cout << "Creating Pores..." << endl;
  createPores();
  cout << "Applying Coordination Number..." << endl;
  applyCoordinationNumber();
  cout << "Cleaning Network..." << endl;
  cleanGenericNetwork();
  cout << "Assigning Radii..." << endl;
  assignRadii();
  cout << "Assigning Lengths..." << endl;
  assignLengths();
  cout << "Distorting Network..." << endl;
  distortNetwork();
  cout << "Assigning Shape Factors..." << endl;
  assignShapeFactors();
  cout << "Assigning Volumes..." << endl;
  assignVolumes();
  cout << "Assigning Conductivities..." << endl;
  assignConductivities();
  cout << "Calculate Network Properties..." << endl;
  setActiveElements();

  if (absolutePermeabilityCalculation) {
    cout << "Absolute Permeabilty Calculation..." << endl;
    solvePressures();
    updateFlows();
    calculatePermeabilityAndPorosity();
  }
}

void network::createNodes() {
  tableOfNodes.resize(Nx);
  for (int i = 0; i < Nx; ++i) {
    tableOfNodes[i].resize(Ny);
    for (int j = 0; j < Ny; ++j) tableOfNodes[i][j].resize(Nz);
  }
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k) {
        tableOfNodes[i][j][k] = new node(i, j, k);
        tableOfAllNodes.push_back(tableOfNodes[i][j][k]);
      }

  totalNodes = Nx * Ny * Nz;

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setId(i + 1);
    if (n->getIndexX() == 0) n->setInlet(true);
    if (n->getIndexX() == Nx - 1) n->setOutlet(true);
    n->setXCoordinate(n->getIndexX() * length);
    n->setYCoordinate(n->getIndexY() * length);
    n->setZCoordinate(n->getIndexZ() * length);
  }

  xEdgeLength = Nx * length;
  yEdgeLength = Ny * length;
  zEdgeLength = Nz * length;
  minNodeRadius = 0;
  maxNodeRadius = 0;
}

void network::createPores() {
  tableOfPoresX.resize(Nx + 1);
  for (int i = 0; i < Nx + 1; ++i) {
    tableOfPoresX[i].resize(Ny + 1);
    for (int j = 0; j < Ny; ++j) tableOfPoresX[i][j].resize(Nz);
  }

  tableOfPoresY.resize(Nx);
  for (int i = 0; i < Nx; ++i) {
    tableOfPoresY[i].resize(Ny + 1);
    for (int j = 0; j < Ny + 1; ++j) tableOfPoresY[i][j].resize(Nz);
  }

  tableOfPoresZ.resize(Nx);
  for (int i = 0; i < Nx; ++i) {
    tableOfPoresZ[i].resize(Ny);
    for (int j = 0; j < Ny; ++j) tableOfPoresZ[i][j].resize(Nz + 1);
  }

  for (int i = 0; i < Nx + 1; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        tableOfPoresX[i][j][k] =
            new pore(getNode(i, j, k), getNode(i - 1, j, k));
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny + 1; ++j)
      for (int k = 0; k < Nz; ++k)
        tableOfPoresY[i][j][k] =
            new pore(getNode(i, j, k), getNode(i, j - 1, k));
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz + 1; ++k)
        tableOfPoresZ[i][j][k] =
            new pore(getNode(i, j, k), getNode(i, j, k - 1));

  totalPores = 3 * Nx * Ny * Nz + Ny * Nz + Nx * Nz + Nx * Ny;

  for (int k = 0; k < Nz; ++k)
    for (int j = 0; j < Ny; ++j)
      for (int i = 0; i < Nx + 1; ++i)
        tableOfAllPores.push_back(tableOfPoresX[i][j][k]);
  for (int k = 0; k < Nz; ++k)
    for (int j = 0; j < Ny + 1; ++j)
      for (int i = 0; i < Nx; ++i)
        tableOfAllPores.push_back(tableOfPoresY[i][j][k]);
  for (int k = 0; k < Nz + 1; ++k)
    for (int j = 0; j < Ny; ++j)
      for (int i = 0; i < Nx; ++i)
        tableOfAllPores.push_back(tableOfPoresZ[i][j][k]);

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
    if (p->getNodeOut() == 0) p->setInlet(true);
    if (p->getNodeIn() == 0) p->setOutlet(true);
  }
  setNeighboors();
}

void network::setNeighboors() {
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    vector<int> connectedPores;
    vector<int> neighboors;
    pore* x = getPoreX(n->getIndexX(), n->getIndexY(), n->getIndexZ());
    pore* xout = getPoreXout(n->getIndexX(), n->getIndexY(), n->getIndexZ());
    pore* y = getPoreY(n->getIndexX(), n->getIndexY(), n->getIndexZ());
    pore* yout = getPoreYout(n->getIndexX(), n->getIndexY(), n->getIndexZ());
    pore* z = getPoreZ(n->getIndexX(), n->getIndexY(), n->getIndexZ());
    pore* zout = getPoreZout(n->getIndexX(), n->getIndexY(), n->getIndexZ());

    connectedPores.push_back(x->getId());
    connectedPores.push_back(xout->getId());
    connectedPores.push_back(y->getId());
    connectedPores.push_back(yout->getId());
    connectedPores.push_back(z->getId());
    connectedPores.push_back(zout->getId());
    n->setConnectedPores(connectedPores);

    neighboors.push_back(x->getNodeOut() == 0 ? 0 : x->getNodeOut()->getId());
    neighboors.push_back(xout->getNodeIn() == 0 ? 0
                                                : xout->getNodeIn()->getId());
    neighboors.push_back(y->getNodeOut() == 0 ? 0 : y->getNodeOut()->getId());
    neighboors.push_back(yout->getNodeIn() == 0 ? 0
                                                : yout->getNodeIn()->getId());
    neighboors.push_back(z->getNodeOut() == 0 ? 0 : z->getNodeOut()->getId());
    neighboors.push_back(zout->getNodeIn() == 0 ? 0
                                                : zout->getNodeIn()->getId());
    n->setNeighboors(neighboors);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    vector<pore*> neighboors;
    if (p->getNodeIn() != 0) {
      const vector<int>& neighboorsIn = p->getNodeIn()->getConnectedPores();
      for (unsigned j = 0; j < neighboorsIn.size(); ++j)
        if (neighboorsIn[j] != p->getId())
          neighboors.push_back(getPore(neighboorsIn[j] - 1));
    }
    if (p->getNodeOut() != 0) {
      const vector<int>& neighboorsOut = p->getNodeOut()->getConnectedPores();
      for (unsigned j = 0; j < neighboorsOut.size(); ++j)
        if (neighboorsOut[j] != p->getId())
          neighboors.push_back(getPore(neighboorsOut[j] - 1));
    }
    p->setNeighboors(neighboors);
  }
}

void network::setNeighboorsForEntireNetwork() {
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    int connectionNumber(0);
    vector<int> neighboorsPores;
    vector<int> neighboorsNodes;

    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      node* nodeIn = p->getNodeIn();
      node* nodeOut = p->getNodeOut();

      if (nodeIn == n || nodeOut == n) {
        neighboorsPores.push_back(p->getId());
        if (nodeIn == n)
          neighboorsNodes.push_back(nodeOut == 0 ? -1 : nodeOut->getId());
        if (nodeOut == n)
          neighboorsNodes.push_back(nodeIn == 0 ? -1 : nodeIn->getId());
        connectionNumber++;
      }
    }
    n->setConnectionNumber(connectionNumber);
    n->setNeighboors(neighboorsNodes);
    n->setConnectedPores(neighboorsPores);
    if (connectionNumber > maxConnectionNumber)
      maxConnectionNumber = connectionNumber;
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    vector<pore*> neighboors;
    if (p->getNodeIn() != 0) {
      const vector<int>& neighboorsIn = p->getNodeIn()->getConnectedPores();
      for (unsigned j = 0; j < neighboorsIn.size(); ++j)
        if (neighboorsIn[j] != p->getId())
          neighboors.push_back(getPore(neighboorsIn[j] - 1));
    }
    if (p->getNodeOut() != 0) {
      const vector<int>& neighboorsOut = p->getNodeOut()->getConnectedPores();
      for (unsigned j = 0; j < neighboorsOut.size(); ++j)
        if (neighboorsOut[j] != p->getId())
          neighboors.push_back(getPore(neighboorsOut[j] - 1));
    }

    p->setNeighboors(neighboors);
  }
}

void network::applyCoordinationNumber() {
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j) {
      pore* zIn = getPoreZ(i, j, 0);
      pore* zOut = getPoreZ(i, j, Nz);
      zIn->setClosed(true);
      zIn->setConductivity(1e-100);
      zIn->getNodeIn()->setConnectionNumber(
          zIn->getNodeIn()->getConnectionNumber() - 1);
      zOut->setClosed(true);
      zOut->setConductivity(1e-100);
      zOut->getNodeOut()->setConnectionNumber(
          zOut->getNodeOut()->getConnectionNumber() - 1);
    }
  for (int i = 0; i < Nx; ++i)
    for (int k = 0; k < Nz; ++k) {
      pore* yIn = getPoreY(i, 0, k);
      pore* yOut = getPoreY(i, Ny, k);
      yIn->setClosed(true);
      yIn->setConductivity(1e-100);
      yIn->getNodeIn()->setConnectionNumber(
          yIn->getNodeIn()->getConnectionNumber() - 1);
      yOut->setClosed(true);
      yOut->setConductivity(1e-100);
      yOut->getNodeOut()->setConnectionNumber(
          yOut->getNodeOut()->getConnectionNumber() - 1);
    }

  if (coordinationNumber < 6 || (coordinationNumber < 4 && Nz == 1)) {
    double totalOpenedPoresSoFar = totalPores - 2 * Nx * Ny - 2 * Nx * Nz;
    int closedPoresNumber =
        (Nz == 1 ? int(totalOpenedPoresSoFar * (1 - coordinationNumber / 4.0))
                 : int(totalPores * (1 - coordinationNumber / 6.0)));
    int i = 0;
    while (i < closedPoresNumber) {
      int index = uniform_int(0, totalPores - 1);
      pore* p = getPore(index);
      if (!p->getClosed()) {
        p->setClosed(true);
        p->setConductivity(1e-100);
        node* in = p->getNodeIn();
        node* out = p->getNodeOut();
        if (in != 0) in->setConnectionNumber(in->getConnectionNumber() - 1);
        if (out != 0) out->setConnectionNumber(out->getConnectionNumber() - 1);
        ++i;
      }
    }
  }

  maxConnectionNumber = (Nz == 1 ? 4 : 6);
}

void network::assignRadii() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      if (radiusDistribution == 1)
        p->setRadius(uniform_real(minRadius, maxRadius));
      if (radiusDistribution == 2)
        p->setRadius(rayleigh(minRadius, maxRadius, rayleighParameter));
      if (radiusDistribution == 3)
        p->setRadius(triangular(minRadius, maxRadius, triangularParameter));
      if (radiusDistribution == 4)
        p->setRadius(weibull(minRadius, maxRadius, normalMuParameter,
                             normalSigmaParameter));
      p->setRadius_sq(pow(p->getRadius(), 2));
    }
  }
}

void network::assignLengths() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setFullLength(length);
      p->setNodeInLength(p->getNodeIn() == 0 ? 0 : p->getNodeIn()->getRadius());
      p->setNodeOutLength(p->getNodeOut() == 0 ? 0
                                               : p->getNodeOut()->getRadius());
      p->setLength(length - p->getNodeInLength() - p->getNodeOutLength() > 0
                       ? length - p->getNodeInLength() - p->getNodeOutLength()
                       : length / 2);
    }
  }
}

void network::distortNetwork() {
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setXCoordinate(n->getXCoordinate() +
                      length * degreeOfDistortion * (-1 + 2 * uniform_real()));
    n->setYCoordinate(n->getYCoordinate() +
                      length * degreeOfDistortion * (-1 + 2 * uniform_real()));
    if (Nz != 1)
      n->setZCoordinate(n->getZCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
  }

  // update pores' lengths
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      if (p->getNodeIn() == 0 || p->getNodeOut() == 0) continue;
      double length = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                   p->getNodeOut()->getXCoordinate(),
                               2) +
                           pow(p->getNodeIn()->getYCoordinate() -
                                   p->getNodeOut()->getYCoordinate(),
                               2) +
                           pow(p->getNodeIn()->getZCoordinate() -
                                   p->getNodeOut()->getZCoordinate(),
                               2));
      p->setFullLength(length);
      p->setLength(length - p->getNodeInLength() - p->getNodeOutLength());
    }
  }
}

void network::assignShapeFactors() {
  // shapeFactor=1./(4*tools::pi());
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

  if (networkSource == 3)
    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      if (!n->getClosed()) n->setShapeFactor(shapeFactor);
    }

  assignShapeFactorConstants();
}

void network::assignShapeFactorConstants() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      if (p->getShapeFactor() <= sqrt(3) / 36.)
        p->setShapeFactorConstant(0.6);
      else if (p->getShapeFactor() <= 1. / 16.)
        p->setShapeFactorConstant(0.5623);
      else
        p->setShapeFactorConstant(0.5);
    }
  }

  if (networkSource == 3)
    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      if (!n->getClosed()) {
        if (n->getShapeFactor() <= sqrt(3) / 36.)
          n->setShapeFactorConstant(0.6);
        else if (n->getShapeFactor() <= 1. / 16.)
          n->setShapeFactorConstant(0.5623);
        else
          n->setShapeFactorConstant(0.5);
      }
    }
}

void network::assignVolumes() {
  totalPoresVolume = 0;
  totalNodesVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = pow(poreVolumeConstant, 2 - poreVolumeExponent) *
                      p->getLength() * pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  if (networkSource == 3)
    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      if (!n->getClosed()) {
        double volume = poreVolumeConstant * n->getLength() *
                        pow(n->getRadius(), poreVolumeExponent) /
                        (4 * n->getShapeFactor()) *
                        pow(10, (6 * poreVolumeExponent - 12));
        n->setVolume(volume);
        totalNodesVolume += volume;
      }
    }
}

void network::assignConductivities() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double conductivity =
          poreConductivityConstant * p->getShapeFactorConstant() *
          pow(p->getRadius(), poreConductivityExponent) /
          (16 * p->getShapeFactor()) / (p->getViscosity() * p->getLength()) *
          pow(10, (6 * poreConductivityExponent - 24));
      ;
      p->setConductivity(conductivity);
    }
  }
}

void network::setActiveElements() {
  totalOpenedPores = 0;
  totalOpenedNodes = 0;

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      totalOpenedNodes++;
    }
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      totalOpenedPores++;
    }
  }
}

void network::assignViscosities() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setViscosity(plasmaViscosity);
  }

  if (networkSource == 3)
    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      n->setViscosity(plasmaViscosity);
    }
}

void network::updateRanking() {
  int rank = 0;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      n->setRank(rank);
      rank++;
    }
  }
}

void network::setNeighboorsForGenericModel() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    vector<pore*> neighboors;
    if (p->getNodeIn() != 0) {
      const vector<int>& neighboorsIn = p->getNodeIn()->getConnectedPores();
      for (unsigned j = 0; j < neighboorsIn.size(); ++j)
        if (neighboorsIn[j] != p->getId())
          neighboors.push_back(getPore(neighboorsIn[j] - 1));
    }
    if (p->getNodeOut() != 0) {
      const vector<int>& neighboorsOut = p->getNodeOut()->getConnectedPores();
      for (unsigned j = 0; j < neighboorsOut.size(); ++j)
        if (neighboorsOut[j] != p->getId())
          neighboors.push_back(getPore(neighboorsOut[j] - 1));
    }
    p->setNeighboors(neighboors);
  }
}

void network::cleanGenericNetwork() {
  // Clean Network from isolated pores

  clusterEverything();
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && !p->getClusterExist()->getSpanning()) {
      node* nodeIn = p->getNodeIn();
      node* nodeOut = p->getNodeOut();
      if (nodeIn != 0) nodeIn->setConnectionNumber(0);
      if (nodeOut != 0) nodeOut->setConnectionNumber(0);
      p->setClosed(true);
    }
  }

  // update coordination number nodes
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    coordinationNumber = 0;
    const vector<int>& neighboors = n->getConnectedPores();
    for (unsigned j = 0; j < neighboors.size(); ++j) {
      pore* p = getPore(neighboors[j] - 1);
      if (!p->getClosed()) coordinationNumber++;
    }
    n->setConnectionNumber(coordinationNumber);
  }

  // update nodes status
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getConnectionNumber() == 0) {
      n->setClosed(true);
    }
    if (n->getConnectionNumber() == 1) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      if (!p->getClosed() && (p->getInlet() || p->getOutlet())) {
        n->setClosed(true);
      }
    }
    if (n->getConnectionNumber() > 1) {
      n->setClosed(false);
    }
  }

  // Update the ranking for the solver
  int rank = 0;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      n->setRank(rank);
      rank++;
    }
  }
}

void network::generateTissue() {
  if (networkSource == 1 || networkSource == 4 || networkSource == 5 ||
      networkSource == 6) {
    meshSizeX = Nx;
    meshSizeY = Ny;
    meshSizeZ = Nz;
  }

  // generate tissue
  double hx = xEdgeLength / meshSizeX;
  double hy = yEdgeLength / meshSizeY;
  double hz = zEdgeLength / meshSizeZ;

  for (int i = 0; i < meshSizeX; ++i)
    for (int j = 0; j < meshSizeY; ++j)
      for (int k = 0; k < meshSizeZ; ++k)
        tableOfAllBlocks.push_back(new block(i, j, k));

  totalBlocks = meshSizeX * meshSizeY * meshSizeZ;

  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    n->setId(i + 1);
    if (n->getX() == 0) n->setInlet(true);
    if (n->getX() == meshSizeX - 1) n->setOutlet(true);
    if (n->getY() == 0) n->setInletY(true);
    if (n->getY() == meshSizeY - 1) n->setOutletY(true);
    if (n->getZ() == 0) n->setInletZ(true);
    if (n->getZ() == meshSizeZ - 1) n->setOutletZ(true);
    n->setXCoordinate(n->getX() * hx);
    n->setYCoordinate(n->getY() * hy);
    n->setZCoordinate(n->getZ() * hz);
    n->setHx(hx);
    n->setHy(hy);
    n->setHz(hz);
    n->setVolume(hx * hy * hz);
    n->setEffectiveVolume(n->getVolume());

    tableOfAllBlocks.push_back(n);
  }

  cout << "Generating Tissue Neighboors..." << endl;

  // Neighbooring blocks
  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    vector<block*> neighboors;
    if (!n->getInlet())
      neighboors.push_back(getBlock(n->getX() - 1, n->getY(), n->getZ()));
    if (!n->getOutlet())
      neighboors.push_back(getBlock(n->getX() + 1, n->getY(), n->getZ()));
    if (!n->getInletY())
      neighboors.push_back(getBlock(n->getX(), n->getY() - 1, n->getZ()));
    if (!n->getOutletY())
      neighboors.push_back(getBlock(n->getX(), n->getY() + 1, n->getZ()));
    if (!n->getInletZ())
      neighboors.push_back(getBlock(n->getX(), n->getY(), n->getZ() - 1));
    if (!n->getOutletZ())
      neighboors.push_back(getBlock(n->getX(), n->getY(), n->getZ() + 1));
    n->setNeighbooringBlocks(neighboors);
  }

  if (networkSource == 1 || networkSource == 4 || networkSource == 5 ||
      networkSource == 6)
    tissueNetworkCollisionAnalysisRegular();
  if (networkSource == 3) tissueNetworkCollisionAnalysis();
}

void network::tissueNetworkCollisionAnalysis() {
  cout << "Tissue Vessels Collision Analysis..." << endl;

  double hx = xEdgeLength / meshSizeX;
  double hy = yEdgeLength / meshSizeY;
  double hz = zEdgeLength / meshSizeZ;
  double coeff = pow(hx, 2) + pow(hy, 2) + pow(hz, 2);
  // Neighbooring vessels
  for (int j = 0; j < totalPores; ++j) {
    pore* p = getPore(j);
    if (!p->getClosed()) {
      set<int> neighbooringBlocks;
      double px2 = p->getNodeIn()->getXCoordinate();
      double py2 = p->getNodeIn()->getYCoordinate();
      double pz2 = p->getNodeIn()->getZCoordinate();
      double px1 = p->getNodeOut()->getXCoordinate();
      double py1 = p->getNodeOut()->getYCoordinate();
      double pz1 = p->getNodeOut()->getZCoordinate();

      double dx = px2 - px1;
      double dy = py2 - py1;
      double dz = pz2 - pz1;

      int cx1 = floor(px1 / hx);
      int cy1 = floor(py1 / hy);
      int cz1 = floor(pz1 / hz);

      int ndx = ceil(abs(dx / hx)) + 1;
      int ndy = ceil(abs(dy / hy)) + 1;
      int ndz = ceil(abs(dz / hz)) + 1;
      int ndrx = ceil(p->getRadius() / hx) + 1;
      int ndry = ceil(p->getRadius() / hy) + 1;
      int ndrz = ceil(p->getRadius() / hz) + 1;

      for (int rx = -ndrx; rx < ndrx; rx++)
        for (int ry = -ndry; ry < ndry; ry++)
          for (int rz = -ndrz; rz < ndrz; rz++)
            for (int x = -abs(ndx); x < abs(ndx); x++)
              for (int y = -abs(ndy); y < abs(ndy); y++)
                for (int z = -abs(ndz); z < abs(ndz); z++) {
                  block* b = getBlock(cx1 + x + rx, cy1 + y + ry, cz1 + z + rz);
                  if (b != 0 && neighbooringBlocks.find(b->getId()) ==
                                    neighbooringBlocks.end()) {
                    int intersection_number(0);
                    for (int xx = 0; xx < 2; xx++)
                      for (int yy = 0; yy < 2; yy++)
                        for (int zz = 0; zz < 2; zz++) {
                          double testX = b->getXCoordinate() +
                                         (xx == 0 ? 1 : -1) * b->getHx() / 2.;
                          double testY = b->getYCoordinate() +
                                         (yy == 0 ? 1 : -1) * b->getHy() / 2.;
                          double testZ = b->getZCoordinate() +
                                         (zz == 0 ? 1 : -1) * b->getHz() / 2.;

                          double pdx, pdy,
                              pdz;  // vector pd from point 1 to test point
                          double dot, dsq;

                          pdx = testX - px1;  // vector from pt1 to test point.
                          pdy = testY - py1;
                          pdz = testZ - pz1;

                          // Dot the d and pd vectors to see if point lies
                          // behind the cylinder cap at pt1.x, pt1.y, pt1.z

                          dot = pdx * dx + pdy * dy + pdz * dz;

                          if (dot > 0.0f && dot < p->getLength_sq()) {
                            // distance squared to the cylinder axis:
                            dsq = (pdx * pdx + pdy * pdy + pdz * pdz) -
                                  dot * dot / p->getLength_sq();

                            if (dsq < p->getRadius_sq()) intersection_number++;
                          }
                        }
                    // check the cube centre
                    {
                      double testX = b->getXCoordinate();
                      double testY = b->getYCoordinate();
                      double testZ = b->getZCoordinate();

                      double pdx, pdy,
                          pdz;  // vector pd from point 1 to test point
                      double dot, dsq;

                      pdx = testX - px1;  // vector from pt1 to test point.
                      pdy = testY - py1;
                      pdz = testZ - pz1;

                      // Dot the d and pd vectors to see if point lies behind
                      // the cylinder cap at pt1.x, pt1.y, pt1.z

                      dot = pdx * dx + pdy * dy + pdz * dz;

                      if (dot > 0.0f && dot < p->getLength_sq()) {
                        // distance squared to the cylinder axis:
                        dsq = (pdx * pdx + pdy * pdy + pdz * pdz) -
                              dot * dot / p->getLength_sq();

                        if (dsq < p->getRadius_sq()) intersection_number++;
                      }
                    }

                    if (intersection_number != 0 && intersection_number != 9) {
                      double area = double(intersection_number) / 9. *
                                    tools::pi() * coeff / 4.;
                      b->neighbooringVesselsArea()[p->getId()] = area;
                      p->neighbooringBlocksArea()[b->getId()] = area;
                    }

                    if (intersection_number == 9) b->setClosed(true);

                    neighbooringBlocks.insert(b->getId());
                  }
                }
    }
  }

  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    if (!n->getClosed()) {
      int intersection_number = n->neighbooringVesselsArea().size();

      if (intersection_number > 0) {
        n->setConnectedToVessel(true);
        n->setVolume((9 - intersection_number) / 9. * n->getVolume());
      }

      if (intersection_number >= 9) n->setClosed(true);
    }
  }

  cout << "Detectings pores not intersecting cubes corners..." << endl;
  // Rescaling the areas

  for (int j = 0; j < totalPores; ++j) {
    pore* p = getPore(j);
    if (!p->getClosed()) {
      double realArea = 2 * tools::pi() * p->getRadius() * p->getLength();

      auto& neigh = p->neighbooringBlocksArea();

      if (neigh.size() == 0) {
        int cx1 = floor(p->getNodeIn()->getXCoordinate() / hx);
        int cy1 = floor(p->getNodeIn()->getYCoordinate() / hy);
        int cz1 = floor(p->getNodeIn()->getZCoordinate() / hz);

        int cx2 = floor(p->getNodeOut()->getXCoordinate() / hx);
        int cy2 = floor(p->getNodeOut()->getYCoordinate() / hy);
        int cz2 = floor(p->getNodeOut()->getZCoordinate() / hz);

        double blocksNumber =
            (abs(cx1 - cx2) + 1) * (abs(cy1 - cy2) + 1) * (abs(cz1 - cz2) + 1);
        for (int x = min(cx1, cx2); x < max(cx1, cx2) + 1; x++)
          for (int y = min(cy1, cy2); y < max(cy1, cy2) + 1; y++)
            for (int z = min(cz1, cz2); z < max(cz1, cz2) + 1; z++) {
              block* b = getBlock(x, y, z);
              if (!b->getClosed()) {
                double area = realArea / blocksNumber;
                b->neighbooringVesselsArea()[p->getId()] = area;
                b->setConnectedToVessel(true);
                p->neighbooringBlocksArea()[b->getId()] = area;
              }
            }
      }
    }
  }

  cout << "Rescaling areas..." << endl;
  // Rescaling the areas

  for (int j = 0; j < totalPores; ++j) {
    pore* p = getPore(j);
    if (!p->getClosed()) {
      double realArea = 2 * tools::pi() * p->getRadius() * p->getLength();
      double area = 0;
      auto& neigh = p->neighbooringBlocksArea();
      for (auto iterator : neigh) {
        int blockID = iterator.first;
        block* bb = getBlock(blockID - 1);
        if (!bb->getClosed()) area += iterator.second;
      }
      if (area == 0) {
        continue;
      }
      double scalingFactor = realArea / area;
      for (auto iterator : neigh) {
        int blockID = iterator.first;
        block* bb = getBlock(blockID - 1);
        if (!bb->getClosed()) {
          double oldArea = iterator.second;
          p->neighbooringBlocksArea()[blockID] = scalingFactor * oldArea;
          bb->neighbooringVesselsArea()[p->getId()] = scalingFactor * oldArea;
        }
      }
    }
  }

  cout << "Rescaling volumes..." << endl;
  // Rescaling the areas

  double volumeOfConnectedBlocks = 0;
  double volumeOfDisconnectedBlocks = 0;
  double totalVesselVolume = 0;
  double voxelVolume = xEdgeLength * yEdgeLength * zEdgeLength;

  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    if (!n->getClosed()) {
      if (n->neighbooringVesselsArea().size() != 0)
        volumeOfConnectedBlocks += n->getVolume();
      else
        volumeOfDisconnectedBlocks += n->getVolume();
    }
  }

  for (int j = 0; j < totalPores; ++j) {
    pore* p = getPore(j);
    if (!p->getClosed()) {
      totalVesselVolume += p->getVolume();
    }
  }

  double adjustingFactor =
      (voxelVolume - volumeOfDisconnectedBlocks - totalVesselVolume) /
      volumeOfConnectedBlocks;

  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    if (!n->getClosed()) {
      if (n->neighbooringVesselsArea().size() != 0) {
        n->setVolume(n->getVolume() * adjustingFactor);
      }
    }
  }
}

void network::tissueNetworkCollisionAnalysisRegular() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      node* nodeIn = p->getNodeIn();
      node* nodeOut = p->getNodeOut();
      if (nodeIn != 0) {
        block* b1 = getBlock(nodeIn->getIndexX(), nodeIn->getIndexY(),
                             nodeIn->getIndexZ());
        b1->neighbooringVesselsArea()[p->getId()] =
            tools::pi() * p->getRadius() * p->getLength();
        b1->setVolume(b1->getVolume() - p->getVolume() / 2);
        b1->setConnectedToVessel(true);
        p->neighbooringBlocksArea()[b1->getId()] =
            tools::pi() * p->getRadius() * p->getLength();
      }
      if (nodeOut != 0) {
        block* b2 = getBlock(nodeOut->getIndexX(), nodeOut->getIndexY(),
                             nodeOut->getIndexZ());
        b2->neighbooringVesselsArea()[p->getId()] =
            tools::pi() * p->getRadius() * p->getLength();
        b2->setVolume(b2->getVolume() - p->getVolume() / 2);
        b2->setConnectedToVessel(true);
        p->neighbooringBlocksArea()[b2->getId()] =
            tools::pi() * p->getRadius() * p->getLength();
      }
    }
  }

  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    if (!n->getClosed()) {
      if (n->getVolume() < 0) n->setClosed(true);
    }
  }
}

void network::setupTissueProperties() {
  if (tissueHomo) {
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        p->setMembranePermeability(PVT);
      }
    }

    for (int i = 0; i < totalBlocks; ++i) {
      block* p = getBlock(i);
      if (!p->getClosed()) {
        p->setDiffusivity(DT);
      }
    }
  }

  // random

  if (tissueRandom) {
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        double dice = uniform_real();
        if (dice < 0.333)
          p->setMembranePermeability(PVT);
        else if (dice < 0.666)
          p->setMembranePermeability(PVT * 10);
        else
          p->setMembranePermeability(PVT * 100.);
      }
    }

    for (int i = 0; i < totalBlocks; ++i) {
      block* p = getBlock(i);
      if (!p->getClosed()) {
        double dice = uniform_real();
        if (dice < 0.333)
          p->setDiffusivity(DT);
        else if (dice < 0.666)
          p->setDiffusivity(DT * 10);
        else
          p->setDiffusivity(DT / 10);
      }
    }
  }

  // circular

  if (tissueCircular) {
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        p->setMembranePermeability(PVT);

        double xc = tissueCircularX * xEdgeLength;
        double yc = tissueCircularY * yEdgeLength;
        double zc = tissueCircularZ * zEdgeLength;

        double x = p->getXCoordinate();
        double y = p->getYCoordinate();
        double z = p->getZCoordinate();

        if (pow(xc - x, 2) + pow(yc - y, 2) + pow(zc - z, 2) <
            pow(tissueCircularR, 2))
          p->setMembranePermeability(PVT * 100);
      }
    }

    for (int i = 0; i < totalBlocks; ++i) {
      block* p = getBlock(i);
      if (!p->getClosed()) {
        p->setDiffusivity(DT);

        double xc = tissueCircularX * xEdgeLength;
        double yc = tissueCircularY * yEdgeLength;
        double zc = tissueCircularZ * zEdgeLength;
        double x = p->getXCoordinate();
        double y = p->getYCoordinate();
        double z = p->getZCoordinate();

        if (pow(xc - x, 2) + pow(yc - y, 2) + pow(zc - z, 2) <
            pow(tissueCircularR, 2))
          p->setDiffusivity(DT / 10);
      }
    }
  }
}
