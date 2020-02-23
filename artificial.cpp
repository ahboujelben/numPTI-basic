/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"

using namespace std;

void network::setupArtificialNetwork() {
  cout << "Constructing Artificial Network.." << endl;
  buildArtificialNetwork();
  cout << "Cleaning Network..." << endl;
  cleanGenericNetwork();
  cout << "Calculate Network Properties..." << endl;
  setActiveElements();

  if (absolutePermeabilityCalculation) {
    cout << "Absolute permeabilty calculation.." << endl;
    solvePressures();
    updateFlows();
    calculatePermeabilityAndPorosity();
  }
}

void network::buildArtificialNetwork() {
  maxConnectionNumber = 0;
  xEdgeLength = 1960. * 1e-6;
  yEdgeLength = 1280. * 1e-6;
  zEdgeLength = 1e-6;
  maxNodeRadius = 0;
  minNodeRadius = 1e-10;
  maxRadius = 20e-6;
  minRadius = 10e-6;
  length = 0;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = 384 + 3200;
  tableOfAllNodes.resize(totalNodes);

  double xcoordinates[] = {200, 170, 120, 80, 50,  25,  10,
                           10,  25,  50,  80, 120, 170, 200};

  tableOfAllNodes[0] = new node(-200e-6, 0, 0);
  tableOfAllNodes[1] = new node(0, 0, 0);

  double abscisse = 0;
  int index = 1;
  int nodeNumber = 2;
  double L = (1710e-6 - 400e-6) / 2;

  for (int i = 0; i < 7; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    index *= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<<nodeNumber<<" "<<abscisse<<" "<<y<<endl;
      nodeNumber++;
    }
  }

  double rabscisse = abscisse;

  for (int i = 0; i < 25; ++i) {
    rabscisse += 10e-6;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(rabscisse, y, 0);
      // cout<<nodeNumber<<" "<<abscisse<<" "<<y<<endl;
      nodeNumber++;
    }
  }

  for (int i = 7; i < 13; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    ;
    rabscisse += xcoordinates[i] * 1e-6;
    ;
    index /= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) +
                 2 * (yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) * j /
                     (index - 1);
      tableOfAllNodes[nodeNumber] = new node(rabscisse, y, 0);
      // cout<<nodeNumber<<" "<<abscisse<<" "<<y<<endl;
      nodeNumber++;
    }
  }

  tableOfAllNodes[nodeNumber] = new node(rabscisse + 200e-6, 0, 0);
  tableOfAllNodes[nodeNumber + 1] = new node(rabscisse + 400e-6, 0, 0);

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setId(i + 1);
    n->setRadius(0);
    n->setYCoordinate(n->getYCoordinate() + yEdgeLength / 2);
    n->setXCoordinate(n->getXCoordinate() + 200e-6);
  }

  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores = 510 + 25 * 128 + 127 * 26;
  tableOfAllPores.resize(totalPores);

  tableOfAllPores[0] = new pore(getNode(1), getNode(0));

  nodeNumber = 1;
  int poreNumber = 1;
  index = 1;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(nodeNumber);

      node* nodeOut1 = getNode(2 * nodeNumber);
      node* nodeOut2 = getNode(2 * nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<poreNumber<<" "<<nodeNumber<<" "<<2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      //   cout<<poreNumber<<" "<<nodeNumber<<" "<<2*nodeNumber+1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  double rnodeNumber = nodeNumber;

  for (int i = 0; i < 25; ++i) {
    for (int j = 0; j < 128; ++j) {
      node* nodeIn = getNode(nodeNumber);
      node* nodeOut = getNode(nodeNumber + 128);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
      // cout<<poreNumber<<" "<<nodeNumber<<" "<<2*nodeNumber<<endl;
      poreNumber++;
      // cout<<nodeNumber<<" "<<abscisse<<" "<<y<<endl;
      nodeNumber++;
    }
  }

  nodeNumber = rnodeNumber;

  for (int i = 0; i < 26; ++i) {
    for (int j = 0; j < 128 - 1; ++j) {
      node* nodeIn = getNode(nodeNumber);
      node* nodeOut = getNode(nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
      // cout<<poreNumber<<" "<<nodeNumber<<" "<<2*nodeNumber<<endl;
      poreNumber++;
      // cout<<nodeNumber<<" "<<abscisse<<" "<<y<<endl;
      nodeNumber++;
    }
    nodeNumber++;
  }

  nodeNumber = 1;
  index = 1;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(totalNodes - 1 - nodeNumber);

      node* nodeOut1 = getNode(totalNodes - 1 - 2 * nodeNumber);
      node* nodeOut2 = getNode(totalNodes - 1 - 2 * nodeNumber - 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      // "<<totalNodes-1-2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      //  "<<totalNodes-1-2*nodeNumber-1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  tableOfAllPores[totalPores - 1] =
      new pore(getNode(totalNodes - 1), getNode(totalNodes - 2));

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(5e-6);
    }
  }
  maxRadius = 6e-6;
  minRadius = 4e-6;

  cout << "Setting neighboors..." << endl;

  // setting neighboors

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
        if (nodeIn == n) neighboorsNodes.push_back(nodeOut->getId());
        if (nodeOut == n) neighboorsNodes.push_back(nodeIn->getId());
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

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setInlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setInlet(true);
      p->setNodeIn(nn);
      p->setNodeOut(n);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }

    if (n->getId() == totalNodes) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setOutlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setOutlet(true);
      p->setNodeIn(n);
      p->setNodeOut(nn);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }
  }

  cout << "Setting Lengths..." << endl;

  degreeOfDistortion = 0.3;
  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1 || n->getId() == 2 || n->getId() == totalNodes ||
        n->getId() == totalNodes - 1)
      continue;
    if (n->getId() < 8 || n->getId() >= totalNodes - 8) {
      n->setXCoordinate(n->getXCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      n->setXCoordinate(n->getXCoordinate() + length / 5 * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length / 5 * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
      length += len;
    }
  }
  length /= totalPores;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(p->getLength() / 200e-6 * 6e-6 + 4e-6);
    }
  }
  maxRadius = 10e-6;
  minRadius = 4e-6;

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  // Assign Conductivities

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

  int closedPoresNumber =
      int(totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (!p->getClosed() && p->getLength() < 12e-6) {
      p->setClosed(true);
      ++i;
    }
  }
}

void network::buildArtificialNetwork2() {
  maxConnectionNumber = 6;
  xEdgeLength = 20000e-6;
  yEdgeLength = 10000e-6;
  zEdgeLength = 1000e-6;
  maxNodeRadius = 100e-6;
  ;
  minNodeRadius = 10e-6;
  maxRadius = 500e-6;
  minRadius = 10e-6;
  length = 1000;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = (5 + 1) * 2 + (10 + 1) * 4 + (50 + 1) * 20 +
                                  (200 + 1) * 80 + (400 + 1) * 160;
  tableOfAllNodes.resize(totalNodes);

  cout << "Layer1..." << endl;
  Ny = 5 + 1;
  Nx = 2;
  length = 2000e-6;
  double oldLength = 4000e-6;
  int nodeNumber(0);

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      n->setId(nodeNumber + 1);
      if (n->getIndexX() == 0) n->setInlet(true);

      n->setXCoordinate(n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  cout << "Layer2..." << endl;

  Ny = 10 + 1;
  Nx = 4;
  length = 1000e-6;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      n->setId(nodeNumber + 1);

      n->setXCoordinate(oldLength + n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  cout << "Layer3..." << endl;

  Ny = 50 + 1;
  Nx = 20;
  length = 200e-6;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      n->setId(nodeNumber + 1);

      n->setXCoordinate(2 * oldLength + n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  cout << "Layer4..." << endl;

  Ny = 200 + 1;
  Nx = 80;
  length = 50e-6;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      n->setId(nodeNumber + 1);

      n->setXCoordinate(3 * oldLength + n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  cout << "Layer5..." << endl;

  Ny = 400 + 1;
  Nx = 160;
  length = 25e-6;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      if (i == Nx - 1) n->setOutlet(true);

      n->setId(nodeNumber + 1);

      n->setXCoordinate(4 * oldLength + n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores = 5 + 1 + (5 + 1) * 2 + (10 + 1) * 4 +
                                  (50 + 1) * 20 + (200 + 1) * 80 +
                                  (400 + 1) * 160 + (5) * 2 + (10) * 4 +
                                  (50) * 20 + (200) * 80 + (400) * 160;
  tableOfAllPores.resize(totalPores);
  double poreNumber(0);
  double nodesSoFar(0);
  nodeNumber = 0;

  cout << "Layer1..." << endl;
  Ny = 5 + 1;
  Nx = 2;
  length = 2000e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == 0)  // inlet pores
      {
        node* nodeOut = 0;
        node* nodeIn = getNode(nodeNumber);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        p->setInlet(true);
        poreNumber++;

        nodeIn->getNeighboors().push_back(-1);
        nodeIn->getConnectedPores().push_back(p->getId());

      } else if (i == Nx - 1) {
        node* nodeOut = getNode(nodeNumber);
        double number = nodesSoFar + j * 2;
        node* nodeIn = getNode(number);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == 0) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }

  cout << "Layer2..." << endl;

  Ny = 10 + 1;
  Nx = 4;
  length = 1000e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(250e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == Nx - 1) {
        node* nodeOut = getNode(nodeNumber);
        double number = nodesSoFar + j * 5;
        node* nodeIn = getNode(number);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(250e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      } else {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(250e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }
  cout << "Layer3..." << endl;

  Ny = 50 + 1;
  Nx = 20;
  length = 200e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(100e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == Nx - 1) {
        node* nodeOut = getNode(nodeNumber);
        double number = nodesSoFar + j * 4;
        node* nodeIn = getNode(number);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(100e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      } else {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(100e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }

  cout << "Layer4..." << endl;

  Ny = 200 + 1;
  Nx = 80;
  length = 50e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(20e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == Nx - 1) {
        node* nodeOut = getNode(nodeNumber);
        double number = nodesSoFar + j * 2;
        node* nodeIn = getNode(number);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(20e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      } else {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(20e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }

  cout << "Layer5..." << endl;

  Ny = 400 + 1;
  Nx = 160;
  length = 25e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(10e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == Nx - 1)  // outlet
      {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = 0;
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(10e-6);
        p->setLength(length);
        p->setOutlet(true);
        poreNumber++;

        nodeOut->getNeighboors().push_back(-1);
        nodeOut->getConnectedPores().push_back(p->getId());
      } else {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(10e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }

  cout << "Setting neighboors..." << endl;

  // setting neighboors

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

  cout << "Create Distortion..." << endl;

  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);

    if (n->getXCoordinate() < 4000e-6) {
      length = 2000e-6;
      if (!n->getInlet()) {
        n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
        n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
      }
    } else if (n->getXCoordinate() < 2 * 4000e-6) {
      length = 1000e-6;
      n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else if (n->getXCoordinate() < 3 * 4000e-6) {
      length = 200e-6;
      n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else if (n->getXCoordinate() < 4 * 4000e-6) {
      length = 50e-6;
      n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      length = 25e-6;
      if (!n->getOutlet()) {
        n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
        n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
      }
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && !p->getInlet() && !p->getOutlet()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
    }
  }

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  cout << "Setting Conductivities..." << endl;

  // Assign Conductivities

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

  cout << "Setting Coordination Number..." << endl;

  int closedPoresNumber =
      max(0.0, totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (p->getXCoordinate() < 4000e-6) continue;
    if (!p->getClosed()) {
      p->setClosed(true);
      ++i;
    }
  }
}
void network::buildArtificialNetwork3()  // branches in 3^n + rectangular part
{
  maxConnectionNumber = 0;
  xEdgeLength = 1960. * 1e-6;
  yEdgeLength = 1280. * 1e-6;
  zEdgeLength = 1e-6;
  maxNodeRadius = 0;
  minNodeRadius = 1e-10;
  maxRadius = 20e-6;
  minRadius = 10e-6;
  length = 0;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = 159 + 4 + 2025;
  tableOfAllNodes.resize(totalNodes);

  double xcoordinates[] = {200, 170, 100, 20, 20, 100, 170, 200};

  tableOfAllNodes[0] = new node(-200e-6, 0, 0);
  tableOfAllNodes[1] = new node(0, 0, 0);

  double abscisse = 0;
  int index = 1;
  int nodeNumber = 2;
  double L = (1710e-6 - 730e-6) / 2;

  for (int i = 0; i < 4; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    index *= 3;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  cout << "index:" << index << endl;
  double rabscisse = abscisse;

  for (int i = 0; i < 25; ++i) {
    rabscisse += 10e-6;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(rabscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  for (int i = 4; i < 7; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    ;
    rabscisse += xcoordinates[i] * 1e-6;
    ;
    index /= 3;

    for (int j = 0; j < index; ++j) {
      double y = -(yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) +
                 2 * (yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) * j /
                     (index - 1);
      tableOfAllNodes[nodeNumber] = new node(rabscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  tableOfAllNodes[nodeNumber] = new node(rabscisse + 200e-6, 0, 0);
  tableOfAllNodes[nodeNumber + 1] = new node(rabscisse + 400e-6, 0, 0);

  cout << "nodeNumber:" << nodeNumber << " "
       << "totalNodes:" << totalNodes << endl;

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setId(i + 1);
    n->setRadius(0);
    n->setYCoordinate(n->getYCoordinate() + yEdgeLength / 2);
    n->setXCoordinate(n->getXCoordinate() + 200e-6);
  }
  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores = 122 + 25 * 81 + 26 * 81 + 14;
  tableOfAllPores.resize(totalPores);

  tableOfAllPores[0] = new pore(getNode(0), getNode(1));

  nodeNumber = 1;
  int poreNumber = 1;
  index = 1;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(nodeNumber);

      node* nodeOut1 = getNode(3 * nodeNumber - 1);
      node* nodeOut2 = getNode(3 * nodeNumber);
      node* nodeOut3 = getNode(3 * nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber+1:"<< 2*nodeNumber+1<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut3);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 3;
  }

  // cout<<"index: " <<index<< endl;

  double rnodeNumber = nodeNumber;

  for (int i = 0; i < 25; ++i) {
    for (int j = 0; j < 81; ++j) {
      node* nodeIn = getNode(nodeNumber);
      node* nodeOut = getNode(nodeNumber + 81);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"nodeNumber+128: "<< nodeNumber+128<<endl;
      poreNumber++;
      nodeNumber++;
    }
  }
  // cout<<"poreNumber after loop:"<< poreNumber<<"nodeNumber after loop:"<<
  // nodeNumber<<endl;

  nodeNumber = rnodeNumber;

  for (int i = 0; i < 25; ++i) {
    for (int j = 0; j < 81 - 1; ++j) {
      node* nodeIn = getNode(nodeNumber);
      node* nodeOut = getNode(nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
      // cout<<"poreNumber:"<<poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*nodeNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      // cout<<"nodeNumber:"<< nodeNumber<<" "<<"abscisse:"<< abscisse<<endl;
      // //<<" "<<"y:"<< y<<endl;
      nodeNumber++;
    }
    nodeNumber++;
  }

  nodeNumber = 1;
  index = 1;
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(totalNodes - 1 - nodeNumber);

      node* nodeOut1 = getNode(totalNodes - 1 - 3 * nodeNumber + 1);
      node* nodeOut2 = getNode(totalNodes - 1 - 3 * nodeNumber);
      node* nodeOut3 = getNode(totalNodes - 1 - 3 * nodeNumber - 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      // "<<totalNodes-1-2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      //  "<<totalNodes-1-2*nodeNumber-1<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut3);
      //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      //  "<<totalNodes-1-2*nodeNumber-1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 3;
  }

  /* index=1;
   for(int i=0;i<1;++i)
   {
       for(int j=0;j<index;++j)
       {
           node* nodeIn=getNode(totalNodes-1-nodeNumber);

           node* nodeOut1=getNode(totalNodes-1-3*nodeNumber+1);
           node* nodeOut2=getNode(totalNodes-1-3*nodeNumber);
           node* nodeOut3=getNode(totalNodes-1-3*nodeNumber-1);

           tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut1);
          // cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
   "<<totalNodes-1-2*nodeNumber<<endl; poreNumber++;
           tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut2);
         //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
   "<<totalNodes-1-2*nodeNumber-1<<endl; poreNumber++;
           tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut3);
         //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
   "<<totalNodes-1-2*nodeNumber-1<<endl; poreNumber++; nodeNumber++;
       }
       index*=3;
   }*/

  /*node* nodeIn=getNode(totalNodes-1-nodeNumber);

  node* nodeOut1=getNode(totalNodes-1-3*nodeNumber+1);
  node* nodeOut2=getNode(totalNodes-1-3*nodeNumber);
  node* nodeOut3=getNode(totalNodes-1-3*nodeNumber-1);

  tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut1);
  tableOfAllPores[poreNumber+1]=new pore(nodeIn,nodeOut2);
  tableOfAllPores[poreNumber+2]=new pore(nodeIn,nodeOut3);*/

  tableOfAllPores[totalPores - 1] =
      new pore(getNode(totalNodes - 2), getNode(totalNodes - 1));

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(5e-4);
    }
  }
  maxRadius = 6e-6;
  minRadius = 4e-6;
  cout << "Setting neighboors..." << endl;

  // setting neighboors

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
        if (nodeIn == n) neighboorsNodes.push_back(nodeOut->getId());
        if (nodeOut == n) neighboorsNodes.push_back(nodeIn->getId());
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

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setInlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setInlet(true);
      p->setNodeIn(nn);
      p->setNodeOut(n);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }

    if (n->getId() == totalNodes) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setOutlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setOutlet(true);
      p->setNodeIn(n);
      p->setNodeOut(nn);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }
  }

  cout << "Setting Lengths..." << endl;

  degreeOfDistortion = 0.3;
  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1 || n->getId() == 2 || n->getId() == totalNodes ||
        n->getId() == totalNodes - 1)
      continue;
    if (n->getId() < 8 || n->getId() >= totalNodes - 8) {
      n->setXCoordinate(n->getXCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      n->setXCoordinate(n->getXCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
      length += len;
    }
  }
  length /= totalPores;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(p->getLength() / 200e-6 * 6e-6 + 4e-6);
    }
  }
  maxRadius = 10e-6;
  minRadius = 4e-6;

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  // Assign Conductivities

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

  int closedPoresNumber =
      int(totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (!p->getClosed() && p->getLength() < 12e-6) {
      p->setClosed(true);
      ++i;
    }
  }
}
void network::buildArtificialNetwork4()  // branches in 2^n
{
  maxConnectionNumber = 0;
  xEdgeLength = 1960. * 1e-6;
  yEdgeLength = 1280. * 1e-6;
  zEdgeLength = 1e-6;
  maxNodeRadius = 0;
  minNodeRadius = 1e-10;
  maxRadius = 20e-6;
  minRadius = 10e-6;
  length = 0;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = 384;
  tableOfAllNodes.resize(totalNodes);

  double xcoordinates[] = {200, 170, 120, 80, 50,  25,  10,
                           10,  25,  50,  80, 120, 170, 10};

  tableOfAllNodes[0] = new node(-200e-6, 0, 0);
  tableOfAllNodes[1] = new node(0, 0, 0);

  double abscisse = 0;
  int index = 1;
  int nodeNumber = 2;
  double L = (1710e-6 - 400e-6) / 2;

  for (int i = 0; i < 7; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    index *= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  cout << "index:" << index << endl;
  /*double rabscisse=abscisse;

  for(int i=0;i<25;++i)
  {
      rabscisse+=10e-6;

      for(int j=0;j<index;++j)
      {
          double
  y=-(abscisse/L)*yEdgeLength/2.+j*(abscisse/L*yEdgeLength)/(index-1);
          tableOfAllNodes[nodeNumber]= new node(rabscisse,y,0);
         //cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
  "<<"y:"<< y<<endl; nodeNumber++;
      }
  }*/

  for (int i = 7; i < 13; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    ;
    // rabscisse+=xcoordinates[i]*1e-6;;
    index /= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) +
                 2 * (yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) * j /
                     (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  tableOfAllNodes[nodeNumber] = new node(abscisse + 200e-6, 0, 0);
  tableOfAllNodes[nodeNumber + 1] = new node(abscisse + 400e-6, 0, 0);

  cout << "nodeNumber:" << nodeNumber << " "
       << "totalNodes:" << totalNodes << endl;

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setId(i + 1);
    n->setRadius(0);
    n->setYCoordinate(n->getYCoordinate() + yEdgeLength / 2);
    n->setXCoordinate(n->getXCoordinate() + 200e-6);
  }

  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores = 510;
  tableOfAllPores.resize(totalPores);

  tableOfAllPores[0] = new pore(getNode(0), getNode(1));

  nodeNumber = 1;
  int poreNumber = 1;
  index = 1;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(nodeNumber);

      node* nodeOut1 = getNode(2 * nodeNumber);
      node* nodeOut2 = getNode(2 * nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber+1:"<< 2*nodeNumber+1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  // cout<<"index: " <<index<< endl;

  /*double rnodeNumber=nodeNumber;

  for(int i=0;i<25;++i)
  {
      for(int j=0;j<128;++j)
      {
          node* nodeIn=getNode(nodeNumber);
          node* nodeOut=getNode(nodeNumber+128);

          tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut);
          //cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
  "<<"nodeNumber+128: "<< nodeNumber+128<<endl; poreNumber++; nodeNumber++;
      }
  }
          //cout<<"poreNumber after loop:"<< poreNumber<<"nodeNumber after
  loop:"<< nodeNumber<<endl;

  nodeNumber=rnodeNumber;

  for(int i=0;i<26;++i)
  {
      for(int j=0;j<128-1;++j)
      {
          node* nodeIn=getNode(nodeNumber);
          node* nodeOut=getNode(nodeNumber+1);

          tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut);
         //cout<<"poreNumber:"<<poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
  "<<"2*nodeNumber:"<< 2*nodeNumber<<endl; poreNumber++;
          //cout<<"nodeNumber:"<< nodeNumber<<" "<<"abscisse:"<< abscisse<<endl;
  //<<" "<<"y:"<< y<<endl; nodeNumber++;
      }
      nodeNumber++;
  }*/

  nodeNumber = 1;
  index = 1;
  for (int i = 0; i < 7; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeOut = getNode(totalNodes - 1 - nodeNumber);

      node* nodeIn1 = getNode(totalNodes - 1 - 2 * nodeNumber);
      node* nodeIn2 = getNode(totalNodes - 1 - 2 * nodeNumber - 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn1, nodeOut);
      // cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      // "<<totalNodes-1-2*nodeNumber<<endl;te
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn2, nodeOut);
      //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      //  "<<totalNodes-1-2*nodeNumber-1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  tableOfAllPores[totalPores - 1] =
      new pore(getNode(totalNodes - 2), getNode(totalNodes - 1));

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(5e-4);
    }
  }
  maxRadius = 6e-6;
  minRadius = 4e-6;

  cout << "Setting neighboors..." << endl;

  // setting neighboors

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
        if (nodeIn == n) neighboorsNodes.push_back(nodeOut->getId());
        if (nodeOut == n) neighboorsNodes.push_back(nodeIn->getId());
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

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setInlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setInlet(true);
      p->setNodeIn(nn);
      p->setNodeOut(n);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }

    if (n->getId() == totalNodes) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setOutlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setOutlet(true);
      p->setNodeIn(n);
      p->setNodeOut(nn);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }
  }

  cout << "Setting Lengths..." << endl;

  degreeOfDistortion = 0.3;
  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1 || n->getId() == 2 || n->getId() == totalNodes ||
        n->getId() == totalNodes - 1)
      continue;
    if (n->getId() < 8 || n->getId() >= totalNodes - 8) {
      n->setXCoordinate(n->getXCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      n->setXCoordinate(n->getXCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
      length += len;
    }
  }
  length /= totalPores;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(p->getLength() / 200e-6 * 6e-6 + 4e-6);
    }
  }
  maxRadius = 10e-6;
  minRadius = 4e-6;

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  // Assign Conductivities

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

  int closedPoresNumber =
      int(totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (!p->getClosed())  // && p->getLength()<12e-6)
    {
      p->setClosed(true);
      ++i;
    }
  }
}
void network::buildArtificialNetwork5()  // small: branches up to 2^2
{
  maxConnectionNumber = 0;
  xEdgeLength = 1960. * 1e-6;
  yEdgeLength = 1280. * 1e-6;
  zEdgeLength = 1e-6;
  maxNodeRadius = 0;
  minNodeRadius = 1e-10;
  maxRadius = 20e-6;
  minRadius = 10e-6;
  length = 0;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = 12;
  tableOfAllNodes.resize(totalNodes);

  double xcoordinates[] = {300, 300, 300, 300};

  tableOfAllNodes[0] = new node(-200e-6, 0, 0);
  tableOfAllNodes[1] = new node(0, 0, 0);

  double abscisse = 0;
  int index = 1;
  int nodeNumber = 2;
  double L = (1710e-6 - 400e-6) / 2;

  for (int i = 0; i < 2; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    index *= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  /* cout<<"index:"<<index<<endl;
   double rabscisse=abscisse;

   for(int i=0;i<25;++i)
   {
       rabscisse+=10e-6;

       for(int j=0;j<index;++j)
       {
           double
   y=-(abscisse/L)*yEdgeLength/2.+j*(abscisse/L*yEdgeLength)/(index-1);
           tableOfAllNodes[nodeNumber]= new node(rabscisse,y,0);
          //cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
   "<<"y:"<< y<<endl; nodeNumber++;
       }
   }*/

  for (int i = 2; i < 3; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    ;
    // rabscisse+=xcoordinates[i]*1e-6;;
    index /= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) +
                 2 * (yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) * j /
                     (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  tableOfAllNodes[nodeNumber] = new node(abscisse + 200e-6, 0, 0);
  tableOfAllNodes[nodeNumber + 1] = new node(abscisse + 400e-6, 0, 0);

  cout << "nodeNumber:" << nodeNumber << " "
       << "totalNodes:" << totalNodes << endl;

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setId(i + 1);
    n->setRadius(0);
    n->setYCoordinate(n->getYCoordinate() + yEdgeLength / 2);
    n->setXCoordinate(n->getXCoordinate() + 200e-6);
  }

  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores = 14;
  tableOfAllPores.resize(totalPores);

  tableOfAllPores[0] = new pore(getNode(0), getNode(1));

  nodeNumber = 1;
  int poreNumber = 1;
  index = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(nodeNumber);

      node* nodeOut1 = getNode(2 * nodeNumber);
      node* nodeOut2 = getNode(2 * nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber+1:"<< 2*nodeNumber+1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  // cout<<"index: " <<index<< endl;

  /*double rnodeNumber=nodeNumber;

  for(int i=0;i<25;++i)
  {
      for(int j=0;j<128;++j)
      {
          node* nodeIn=getNode(nodeNumber);
          node* nodeOut=getNode(nodeNumber+128);

          tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut);
          //cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
  "<<"nodeNumber+128: "<< nodeNumber+128<<endl; poreNumber++; nodeNumber++;
      }
  }
          //cout<<"poreNumber after loop:"<< poreNumber<<"nodeNumber after
  loop:"<< nodeNumber<<endl;

  nodeNumber=rnodeNumber;*/

  /*for(int i=0;i<2;++i)
  {
      for(int j=0;j<4-1;++j)
      {
          node* nodeIn=getNode(nodeNumber);
          node* nodeOut=getNode(nodeNumber+1);

          tableOfAllPores[poreNumber]=new pore(nodeIn,nodeOut);
         //cout<<"poreNumber:"<<poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
  "<<"2*nodeNumber:"<< 2*nodeNumber<<endl; poreNumber++;
          //cout<<"nodeNumber:"<< nodeNumber<<" "<<"abscisse:"<< abscisse<<endl;
  //<<" "<<"y:"<< y<<endl; nodeNumber++;
      }
      nodeNumber++;
  }
*/

  nodeNumber = 1;
  index = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeOut = getNode(totalNodes - 1 - nodeNumber);

      node* nodeIn1 = getNode(totalNodes - 1 - 2 * nodeNumber);
      node* nodeIn2 = getNode(totalNodes - 1 - 2 * nodeNumber - 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn1, nodeOut);
      // cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      // "<<totalNodes-1-2*nodeNumber<<endl;te
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn2, nodeOut);
      //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      //  "<<totalNodes-1-2*nodeNumber-1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  tableOfAllPores[totalPores - 1] =
      new pore(getNode(totalNodes - 2), getNode(totalNodes - 1));

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(5e-4);
    }
  }
  maxRadius = 6e-6;
  minRadius = 4e-6;

  cout << "Setting neighboors..." << endl;

  // setting neighboors

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
        if (nodeIn == n) neighboorsNodes.push_back(nodeOut->getId());
        if (nodeOut == n) neighboorsNodes.push_back(nodeIn->getId());
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

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setInlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setInlet(true);
      p->setNodeIn(nn);
      p->setNodeOut(n);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }

    if (n->getId() == totalNodes) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setOutlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setOutlet(true);
      p->setNodeIn(n);
      p->setNodeOut(nn);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }
  }

  cout << "Setting Lengths..." << endl;

  degreeOfDistortion = 0.3;
  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1 || n->getId() == 2 || n->getId() == totalNodes ||
        n->getId() == totalNodes - 1)
      continue;
    if (n->getId() < 8 || n->getId() >= totalNodes - 8) {
      n->setXCoordinate(n->getXCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      n->setXCoordinate(n->getXCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
      length += len;
    }
  }
  length /= totalPores;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(p->getLength() / 200e-6 * 6e-6 + 4e-6);
    }
  }
  maxRadius = 10e-6;
  minRadius = 4e-6;

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  // Assign Conductivities

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

  int closedPoresNumber =
      int(totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (!p->getClosed())  //&& p->getLength()<12e-6)
    {
      p->setClosed(true);
      ++i;
    }
  }
}
void network::buildArtificialNetwork6()  // small: branches up to 2^2 +
                                         // rectangular part
{
  maxConnectionNumber = 0;
  xEdgeLength = 1960. * 1e-6;
  yEdgeLength = 1280. * 1e-6;
  zEdgeLength = 1e-6;
  maxNodeRadius = 0;
  minNodeRadius = 1e-10;
  maxRadius = 20e-6;
  minRadius = 10e-6;
  length = 0;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = 12 + 2 * 4;
  tableOfAllNodes.resize(totalNodes);

  double xcoordinates[] = {300, 300, 300, 300};

  tableOfAllNodes[0] = new node(-200e-6, 0, 0);
  tableOfAllNodes[1] = new node(0, 0, 0);

  double abscisse = 0;
  int index = 1;
  int nodeNumber = 2;
  double L = (1710e-6 - 400e-6) / 2;

  for (int i = 0; i < 2; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    index *= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(abscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  cout << "index:" << index << endl;
  double rabscisse = abscisse;

  for (int i = 0; i < 2; ++i) {
    rabscisse += 100e-6;

    for (int j = 0; j < index; ++j) {
      double y = -(abscisse / L) * yEdgeLength / 2. +
                 j * (abscisse / L * yEdgeLength) / (index - 1);
      tableOfAllNodes[nodeNumber] = new node(rabscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  for (int i = 2; i < 3; ++i) {
    abscisse += xcoordinates[i] * 1e-6;
    ;
    rabscisse += xcoordinates[i] * 1e-6;
    ;
    index /= 2;

    for (int j = 0; j < index; ++j) {
      double y = -(yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) +
                 2 * (yEdgeLength - (abscisse / (2 * L)) * yEdgeLength) * j /
                     (index - 1);
      tableOfAllNodes[nodeNumber] = new node(rabscisse, y, 0);
      // cout<< "nodeNumber:"<<nodeNumber<<" "<<"abscisse:"<< abscisse<<"
      // "<<"y:"<< y<<endl;
      nodeNumber++;
    }
  }

  tableOfAllNodes[nodeNumber] = new node(rabscisse + 200e-6, 0, 0);
  tableOfAllNodes[nodeNumber + 1] = new node(rabscisse + 400e-6, 0, 0);

  cout << "nodeNumber:" << nodeNumber << " "
       << "totalNodes:" << totalNodes << endl;

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    n->setId(i + 1);
    n->setRadius(0);
    n->setYCoordinate(n->getYCoordinate() + yEdgeLength / 2);
    n->setXCoordinate(n->getXCoordinate() + 200e-6);
  }

  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores = 14 + 2 * 4 + 2 * 3;
  tableOfAllPores.resize(totalPores);

  tableOfAllPores[0] = new pore(getNode(0), getNode(1));

  nodeNumber = 1;
  int poreNumber = 1;
  index = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeIn = getNode(nodeNumber);

      node* nodeOut1 = getNode(2 * nodeNumber);
      node* nodeOut2 = getNode(2 * nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut1);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut2);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*noedNumber+1:"<< 2*nodeNumber+1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  // cout<<"index: " <<index<< endl;

  double rnodeNumber = nodeNumber;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4; ++j) {
      node* nodeIn = getNode(nodeNumber);
      node* nodeOut = getNode(nodeNumber + 4);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
      // cout<<"poreNumber:"<< poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"nodeNumber+128: "<< nodeNumber+128<<endl;
      poreNumber++;
      nodeNumber++;
    }
  }
  // cout<<"poreNumber after loop:"<< poreNumber<<"nodeNumber after loop:"<<
  // nodeNumber<<endl;

  nodeNumber = rnodeNumber;

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < 4 - 1; ++j) {
      node* nodeIn = getNode(nodeNumber);
      node* nodeOut = getNode(nodeNumber + 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
      // cout<<"poreNumber:"<<poreNumber<<" "<<"nodeNumber:"<< nodeNumber<<"
      // "<<"2*nodeNumber:"<< 2*nodeNumber<<endl;
      poreNumber++;
      // cout<<"nodeNumber:"<< nodeNumber<<" "<<"abscisse:"<< abscisse<<endl;
      // //<<" "<<"y:"<< y<<endl;
      nodeNumber++;
    }
    nodeNumber++;
  }

  nodeNumber = 1;
  index = 1;
  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < index; ++j) {
      node* nodeOut = getNode(totalNodes - 1 - nodeNumber);

      node* nodeIn1 = getNode(totalNodes - 1 - 2 * nodeNumber);
      node* nodeIn2 = getNode(totalNodes - 1 - 2 * nodeNumber - 1);

      tableOfAllPores[poreNumber] = new pore(nodeIn1, nodeOut);
      // cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      // "<<totalNodes-1-2*nodeNumber<<endl;te
      poreNumber++;
      tableOfAllPores[poreNumber] = new pore(nodeIn2, nodeOut);
      //  cout<<poreNumber<<" "<<totalNodes-1-nodeNumber<<"
      //  "<<totalNodes-1-2*nodeNumber-1<<endl;
      poreNumber++;
      nodeNumber++;
    }
    index *= 2;
  }

  tableOfAllPores[totalPores - 1] =
      new pore(getNode(totalNodes - 2), getNode(totalNodes - 1));

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(5e-4);
    }
  }
  maxRadius = 6e-6;
  minRadius = 4e-6;

  cout << "Setting neighboors..." << endl;

  // setting neighboors

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
        if (nodeIn == n) neighboorsNodes.push_back(nodeOut->getId());
        if (nodeOut == n) neighboorsNodes.push_back(nodeIn->getId());
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

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setInlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setInlet(true);
      p->setNodeIn(nn);
      p->setNodeOut(n);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }

    if (n->getId() == totalNodes) {
      pore* p = getPore(n->getConnectedPores()[0] - 1);
      p->setOutlet(true);
      node* nn = getNode(n->getNeighboors()[0] - 1);
      nn->setOutlet(true);
      p->setNodeIn(n);
      p->setNodeOut(nn);
      n->setClosed(true);
      n->setConductivity(1e-100);
    }
  }

  cout << "Setting Lengths..." << endl;

  degreeOfDistortion = 0.3;
  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (n->getId() == 1 || n->getId() == 2 || n->getId() == totalNodes ||
        n->getId() == totalNodes - 1)
      continue;
    if (n->getId() < 8 || n->getId() >= totalNodes - 8) {
      n->setXCoordinate(n->getXCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + 5 * length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      n->setXCoordinate(n->getXCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() +
                        length / 5 * 0.42 * (-1 + 2 * uniform_real()));
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
      length += len;
    }
  }
  length /= totalPores;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setRadius(p->getLength() / 200e-6 * 6e-6 + 4e-6);
    }
  }
  maxRadius = 10e-6;
  minRadius = 4e-6;

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  // Assign Conductivities

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

  int closedPoresNumber =
      int(totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (!p->getClosed())  //&& p->getLength()<12e-6)
    {
      p->setClosed(true);
      ++i;
    }
  }
}
void network::buildArtificialNetwork7()  // small: 2 layers of different
                                         // densities
{
  maxConnectionNumber = 6;
  xEdgeLength = 11000e-6;
  yEdgeLength = 15500e-6;
  zEdgeLength = 1000e-6;
  maxNodeRadius = 100e-6;
  ;
  minNodeRadius = 10e-6;
  maxRadius = 500e-6;
  minRadius = 10e-6;
  length = 1000;
  Nz = 1;

  cout << "Loading nodes..." << endl;

  // Creating Nodes
  totalOpenedNodes = totalNodes = (5 + 1) * 2 + (10 + 1) * 4;
  tableOfAllNodes.resize(totalNodes);

  cout << "Layer1..." << endl;
  Ny = 5 + 1;
  Nx = 2;
  length = 2000e-6;
  double oldLength = 4000e-6;
  int nodeNumber(0);

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      n->setId(nodeNumber + 1);
      if (n->getIndexX() == 0) n->setInlet(true);

      n->setXCoordinate(n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length + 2.5e-3);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  cout << "Layer2..." << endl;

  Ny = 10 + 1;
  Nx = 4;
  length = 1000e-6;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      tableOfAllNodes[nodeNumber] = new node(i, j, 0);
      node* n = tableOfAllNodes[nodeNumber];

      if (i == Nx - 1) n->setOutlet(true);
      n->setId(nodeNumber + 1);

      n->setXCoordinate(oldLength + n->getIndexX() * length);
      n->setYCoordinate(n->getIndexY() * length + 2.5e-3);
      n->setZCoordinate(n->getIndexZ() * length);

      n->setRadius(minNodeRadius);

      nodeNumber++;
    }
  }

  ///////////

  cout << "Loading pores..." << endl;

  // Creating Pores

  totalOpenedPores = totalPores =
      5 + 1 + (5 + 1) * 2 + (10 + 1) * 4 + (5) * 2 + (10) * 4;
  tableOfAllPores.resize(totalPores);
  double poreNumber(0);
  double nodesSoFar(0);
  nodeNumber = 0;

  cout << "Layer1..." << endl;
  Ny = 5 + 1;
  Nx = 2;
  length = 2000e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == 0)  // inlet pores
      {
        node* nodeOut = 0;
        node* nodeIn = getNode(nodeNumber);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        p->setInlet(true);
        poreNumber++;

        nodeIn->getNeighboors().push_back(-1);
        nodeIn->getConnectedPores().push_back(p->getId());

      } else if (i == Nx - 1) {
        node* nodeOut = getNode(nodeNumber);
        double number = nodesSoFar + j * 2;
        node* nodeIn = getNode(number);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == 0) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(500e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }

  cout << "Layer2..." << endl;

  Ny = 10 + 1;
  Nx = 4;
  length = 1000e-6;
  nodesSoFar += Nx * Ny;

  for (int i = 0; i < Nx; ++i) {
    for (int j = 0; j < Ny; ++j) {
      if (j != Ny - 1) {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + 1);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(250e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      if (i == Nx - 1)  // Outlet
      {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = 0;
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(250e-6);
        p->setLength(length);
        p->setOutlet(true);
        poreNumber++;

        nodeOut->getNeighboors().push_back(-1);
        nodeOut->getConnectedPores().push_back(p->getId());
      } else {
        node* nodeOut = getNode(nodeNumber);
        node* nodeIn = getNode(nodeNumber + Ny);
        tableOfAllPores[poreNumber] = new pore(nodeIn, nodeOut);
        pore* p = tableOfAllPores[poreNumber];
        p->setId(poreNumber + 1);
        p->setRadius(250e-6);
        p->setLength(length);
        poreNumber++;

        nodeIn->getNeighboors().push_back(nodeOut->getId());
        nodeIn->getConnectedPores().push_back(p->getId());
        nodeOut->getNeighboors().push_back(nodeIn->getId());
        nodeOut->getConnectedPores().push_back(p->getId());
      }

      nodeNumber++;
    }
  }

  cout << "Setting neighboors..." << endl;

  // setting neighboors

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

  cout << "Create Distortion..." << endl;

  length = 50e-6;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);

    if (n->getXCoordinate() < 4000e-6) {
      length = 2000e-6;
      if (!n->getInlet()) {
        n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
        n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
      }
    } else if (n->getXCoordinate() < 2 * 4000e-6) {
      length = 1000e-6;
      n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else if (n->getXCoordinate() < 3 * 4000e-6) {
      length = 200e-6;
      n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else if (n->getXCoordinate() < 4 * 4000e-6) {
      length = 50e-6;
      n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
      n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                  (-1 + 2 * uniform_real()));
    } else {
      length = 25e-6;
      if (!n->getOutlet()) {
        n->setXCoordinate(n->getXCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
        n->setYCoordinate(n->getYCoordinate() + length * degreeOfDistortion *
                                                    (-1 + 2 * uniform_real()));
      }
    }
  }

  // Assign Length
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && !p->getInlet() && !p->getOutlet()) {
      double len = sqrt(pow(p->getNodeIn()->getXCoordinate() -
                                p->getNodeOut()->getXCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getYCoordinate() -
                                p->getNodeOut()->getYCoordinate(),
                            2) +
                        pow(p->getNodeIn()->getZCoordinate() -
                                p->getNodeOut()->getZCoordinate(),
                            2));
      p->setLength(len);
      p->setFullLength(len);
    }
  }

  cout << "Setting Shape Factor..." << endl;

  // Assign shape factor

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) p->setShapeFactor(shapeFactor);
  }

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

  cout << "Setting Volumes..." << endl;

  // Assign Volumes

  totalPoresVolume = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      double volume = poreVolumeConstant * p->getLength() *
                      pow(p->getRadius(), poreVolumeExponent) /
                      (4 * p->getShapeFactor()) *
                      pow(10, (6 * poreVolumeExponent - 12));
      p->setVolume(volume);
      totalPoresVolume += volume;
    }
  }

  cout << "Setting Conductivities..." << endl;

  // Assign Conductivities

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

  cout << "Setting Coordination Number..." << endl;

  int closedPoresNumber =
      max(0.0, totalOpenedPores * (1 - coordinationNumber / 4.0));
  int i = 0;
  while (i < closedPoresNumber) {
    int index = uniform_int(0, totalPores - 1);
    pore* p = getPore(index);
    if (p->getXCoordinate() < 4000e-6) continue;
    if (!p->getClosed()) {
      p->setClosed(true);
      ++i;
    }
  }
}
