/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

using namespace std;

void network::createParentVessel() {
  length = 10e-6;

  xEdgeLength = Nx * length;
  yEdgeLength = Ny * length;
  zEdgeLength = Nz * length;
  h_sq = pow(1 / double(Nx), 2);

  cout << "Creating Parent Vessel..." << endl;
  totalPores = totalOpenedPores = Ny + 1;
  totalNodes = totalOpenedNodes = Ny;

  maxConnectionNumber = (Nz == 1 ? 4 : 6);

  cout << "Creating Nodes..." << endl;

  tableOfNodes.resize(Nx);
  for (int i = 0; i < Nx; ++i) {
    tableOfNodes[i].resize(Ny);
    for (int j = 0; j < Ny; ++j) tableOfNodes[i][j].resize(Nz);
  }
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k) {
        tableOfNodes[i][j][k] = 0;
      }

  for (int i = 0; i < Ny; i++) {
    node* n = new node(0, i, 0);
    tableOfAllNodes.push_back(n);
    tableOfNodes[0][i][0] = n;
    n->setId(i + 1);
    n->setXCoordinate(n->getIndexX() * length);
    n->setYCoordinate(n->getIndexY() * length);
    n->setZCoordinate(n->getIndexZ() * length);
  }

  cout << "Creating Pores..." << endl;
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
      for (int k = 0; k < Nz; ++k) tableOfPoresX[i][j][k] = 0;
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny + 1; ++j)
      for (int k = 0; k < Nz; ++k) tableOfPoresY[i][j][k] = 0;
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz + 1; ++k) tableOfPoresZ[i][j][k] = 0;

  pore* inletPore = new pore(getNode(0), 0);
  tableOfAllPores.push_back(inletPore);
  tableOfPoresX[0][0][0] = inletPore;
  pore* outletPore = new pore(0, getNode(Ny - 1));
  tableOfAllPores.push_back(outletPore);
  tableOfPoresX[0][Ny - 1][0] = outletPore;

  for (int i = 0; i < Ny - 1; i++) {
    pore* p = new pore(getNode(i + 1), getNode(i));
    tableOfAllPores.push_back(p);
    tableOfPoresY[0][i + 1][0] = p;
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
    p->setRadius(14e-6);
    if (p->getNodeOut() == 0) p->setInlet(true);
    if (p->getNodeIn() == 0) p->setOutlet(true);
    p->setParentVessel(true);
    p->setHDConcentration(1);
  }

  cout << "Setting Neighboors..." << endl;
  setNeighboorsForEntireNetwork();
  cout << "Cleaning Network..." << endl;
  cleanGenericNetwork();
  cout << "Assigning Lengths..." << endl;
  assignLengths();
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

void network::createParentVessel3D() {
  length = 10e-6;

  xEdgeLength = Nx * length;
  yEdgeLength = Ny * length;
  zEdgeLength = Nz * length;
  h_sq = pow(1 / double(Nx), 2);

  cout << "Creating Parent Vessel..." << endl;
  totalPores = totalOpenedPores = (Ny + 1) * 2 + 2;
  totalNodes = totalOpenedNodes = Ny * 2;

  maxConnectionNumber = (Nz == 1 ? 4 : 6);

  cout << "Creating Nodes..." << endl;

  tableOfNodes.resize(Nx);
  for (int i = 0; i < Nx; ++i) {
    tableOfNodes[i].resize(Ny);
    for (int j = 0; j < Ny; ++j) tableOfNodes[i][j].resize(Nz);
  }
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k) {
        tableOfNodes[i][j][k] = 0;
      }

  for (int i = 0; i < Ny; i++) {
    node* n = new node(0, i, Nz / 3);
    tableOfAllNodes.push_back(n);
    tableOfNodes[0][i][Nz / 3] = n;
    n->setId(i + 1);
    n->setXCoordinate(n->getIndexX() * length);
    n->setYCoordinate(n->getIndexY() * length);
    n->setZCoordinate(n->getIndexZ() * length);
  }

  for (int i = 0; i < Ny; i++) {
    node* n = new node(0, i, 2 * Nz / 3);
    tableOfAllNodes.push_back(n);
    tableOfNodes[0][i][2 * Nz / 3] = n;
    n->setId(i + 1 + Ny);
    n->setXCoordinate(n->getIndexX() * length);
    n->setYCoordinate(n->getIndexY() * length);
    n->setZCoordinate(n->getIndexZ() * length);
  }

  cout << "Creating Pores..." << endl;
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
      for (int k = 0; k < Nz; ++k) tableOfPoresX[i][j][k] = 0;
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny + 1; ++j)
      for (int k = 0; k < Nz; ++k) tableOfPoresY[i][j][k] = 0;
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz + 1; ++k) tableOfPoresZ[i][j][k] = 0;

  pore* inletPore = new pore(getNode(0), 0);
  tableOfAllPores.push_back(inletPore);
  tableOfPoresX[0][0][Nz / 3] = inletPore;
  pore* outletPore = new pore(0, getNode(Ny - 1));
  tableOfAllPores.push_back(outletPore);
  tableOfPoresX[0][Ny - 1][Nz / 3] = outletPore;

  inletPore = new pore(getNode(0 + Ny), 0);
  tableOfAllPores.push_back(inletPore);
  tableOfPoresX[0][0][2 * Nz / 3] = inletPore;
  outletPore = new pore(0, getNode(Ny - 1 + Ny));
  tableOfAllPores.push_back(outletPore);
  tableOfPoresX[0][Ny - 1][2 * Nz / 3] = outletPore;

  pore* xPore = new pore(getNode(0), getNode(Ny));
  xPore->setVesselType(1);
  tableOfAllPores.push_back(xPore);

  pore* xPore2 = new pore(getNode(Ny - 1), getNode(Ny - 1 + Ny));
  tableOfAllPores.push_back(xPore2);
  xPore2->setVesselType(1);

  for (int i = 0; i < Ny - 1; i++) {
    pore* p = new pore(getNode(i + 1), getNode(i));
    tableOfAllPores.push_back(p);
    tableOfPoresY[0][i + 1][Nz / 3] = p;
  }

  for (int i = 0; i < Ny - 1; i++) {
    pore* p = new pore(getNode(i + 1 + Ny), getNode(i + Ny));
    tableOfAllPores.push_back(p);
    tableOfPoresY[0][i + 1][2 * Nz / 3] = p;
  }

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    p->setId(i + 1);
    p->setRadius(14e-6);
    if (p->getNodeOut() == 0) p->setInlet(true);
    if (p->getNodeIn() == 0) p->setOutlet(true);
    p->setParentVessel(true);
    p->setHDConcentration(1);
  }

  cout << "Setting Neighboors..." << endl;
  setNeighboorsForEntireNetwork();
  cout << "Cleaning Network..." << endl;
  cleanGenericNetwork();
  cout << "Assigning Lengths..." << endl;
  assignLengths();
  xPore->setLength(length * Nz / 3);
  xPore2->setLength(length * Nz / 3);
  xPore->setFullLength(length * Nz / 3);
  xPore2->setFullLength(length * Nz / 3);
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

void network::createRetinaParentVessels() {
  length = 20e-6;

  Nx = 220;
  Ny = 220;
  Nz = 1;

  xEdgeLength = Nx * length;
  yEdgeLength = Ny * length;
  zEdgeLength = Nz * length;
  h_sq = pow(1 / double(Nx), 2);

  cout << "Creating Parent Vessel..." << endl;
  totalPores = totalOpenedPores = 10 + 2;
  totalNodes = totalOpenedNodes = 10 + 2;

  maxConnectionNumber = 4;

  // Initialising tables

  tableOfNodes.resize(Nx);
  for (int i = 0; i < Nx; ++i) {
    tableOfNodes[i].resize(Ny);
    for (int j = 0; j < Ny; ++j) tableOfNodes[i][j].resize(Nz);
  }
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k) {
        tableOfNodes[i][j][k] = 0;
      }

  tableOfPoresX.resize(Nx + 1);
  for (int i = 0; i < Nx + 1; ++i) {
    tableOfPoresX[i].resize(Ny);
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
      for (int k = 0; k < Nz; ++k) tableOfPoresX[i][j][k] = 0;

  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny + 1; ++j)
      for (int k = 0; k < Nz; ++k) tableOfPoresY[i][j][k] = 0;
  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz + 1; ++k) tableOfPoresZ[i][j][k] = 0;

  // Creating elements

  cout << "Creating Central nodes..." << endl;
  int cx(Nx / 2), cy(Ny / 2);
  node* ncenter1 = new node(cx, cy, 0);
  tableOfAllNodes.push_back(ncenter1);
  ncenter1->setId(1);
  ncenter1->setXCoordinate(ncenter1->getIndexX() * length);
  ncenter1->setYCoordinate(ncenter1->getIndexY() * length);
  ncenter1->setZCoordinate(ncenter1->getIndexZ() * length);

  node* ncenter2 = new node(cx, cy, 0);
  tableOfAllNodes.push_back(ncenter2);
  ncenter2->setId(2);
  ncenter2->setXCoordinate(ncenter2->getIndexX() * length);
  ncenter2->setYCoordinate(ncenter2->getIndexY() * length);
  ncenter2->setZCoordinate(ncenter2->getIndexZ() * length);

  cout << "Creating inlet/outlet vessels..." << endl;
  // inlet nad outlet pore
  pore* inletPore = new pore(ncenter1, 0);
  tableOfAllPores.push_back(inletPore);
  inletPore->setId(1);
  inletPore->setRadius(14e-6);
  inletPore->setParentVessel(true);
  inletPore->setHDConcentration(1);
  inletPore->setLength(length);
  inletPore->setFullLength(length);
  inletPore->setInlet(true);

  pore* outletPore = new pore(0, ncenter2);
  tableOfAllPores.push_back(outletPore);
  outletPore->setId(2);
  outletPore->setRadius(14e-6);
  outletPore->setParentVessel(true);
  outletPore->setHDConcentration(1);
  outletPore->setLength(length);
  outletPore->setFullLength(length);
  outletPore->setOutlet(true);

  cout << "Creating feeding vessels..." << endl;
  int nodalRadius = 15;
  int currentNodeIndex = 2;
  int currentPoreIndex = 2;
  for (int i = 0; i < 10; ++i) {
    double theta = i * tools::pi() / 5;
    double ix = int(nodalRadius * cos(theta));
    double iy = int(nodalRadius * sin(theta));

    // create a node at the bounday of the circle
    node* n = new node(Nx / 2 + ix, Ny / 2 + iy, 0);
    tableOfAllNodes.push_back(n);
    tableOfNodes[Nx / 2 + ix][Ny / 2 + iy][0] = n;
    n->setId(++currentNodeIndex);
    n->setXCoordinate(n->getIndexX() * length);
    n->setYCoordinate(n->getIndexY() * length);
    n->setZCoordinate(n->getIndexZ() * length);

    // create the associated pore
    node* ncenter = i % 2 ? ncenter1 : ncenter2;
    pore* p = new pore(n, ncenter);
    tableOfAllPores.push_back(p);
    p->setId(++currentPoreIndex);
    p->setRadius(10e-6);
    p->setParentVessel(true);
    p->setHDConcentration(1);
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

  cout << "Setting Neighboors..." << endl;
  setNeighboorsForEntireNetwork();
  cout << "Assigning Lengths..." << endl;
  assignLengths();
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
