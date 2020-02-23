/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

using namespace std;

network::network(QObject *parent) : QObject(parent) { reset(); }

network::~network() {
  for (int i = 0; i < totalPores; ++i) delete tableOfAllPores[i];

  for (int i = 0; i < totalNodes; ++i) delete tableOfAllNodes[i];

  for (int i = 0; i < totalBlocks; ++i) delete tableOfAllBlocks[i];

  for (int i = 0; i < totalParticles; ++i) delete tableOfParticles[i];

  if (!existClusters.empty())
    for (unsigned i = 0; i < existClusters.size(); ++i) delete existClusters[i];

  tools::cleanVideosFolder();
}

void network::destroy() {
  for (int i = 0; i < totalPores; ++i) delete tableOfAllPores[i];

  for (int i = 0; i < totalNodes; ++i) delete tableOfAllNodes[i];

  for (int i = 0; i < totalBlocks; ++i) delete tableOfAllBlocks[i];

  for (int i = 0; i < totalParticles; ++i) delete tableOfParticles[i];

  totalPores = 0;
  totalNodes = 0;
  totalBlocks = 0;
  totalParticles = 0;
  tableOfAllPores.clear();
  tableOfAllNodes.clear();
  tableOfAllBlocks.clear();
  tableOfNodes.clear();
  tableOfPoresX.clear();
  tableOfPoresY.clear();
  tableOfPoresZ.clear();
  tableOfParticles.clear();

  if (!existClusters.empty())
    for (unsigned i = 0; i < existClusters.size(); ++i) delete existClusters[i];

  existClusters.clear();

  tools::cleanVideosFolder();
}

void network::reset() {
  pressureIn = 1;
  pressureOut = 0;

  totalPores = 0;
  totalNodes = 0;
  totalBlocks = 0;
  totalParticles = 0;

  record = false;
  videoRecording = false;
  ready = false;
  cancel = false;
  simulationRunning = false;

  shapeFactor = 1 / (4 * tools::pi());
}

void network::setupModel() {
  if (ready) {
    ready = false;
    destroy();
  }

  reset();
  loadNetworkData();

  gen.seed(seed);

  if (networkSource == 1)  // regular
    setupRegularModel();

  if (networkSource == 2)  // artificial
    setupArtificialNetwork();

  if (networkSource == 4)  // parent vessel tumour
  {
    if (Nz < 5)
      createParentVessel();
    else
      createParentVessel3D();
  }

  if (networkSource == 5)  // parent vessel retina
    createRetinaParentVessels();

  // Build tissue
  if (buildTissue) generateTissue();

  // Notify display to update graphics
  ready = true;
  emitPlotSignal();
}

void network::runSimulation() {
  tools::cleanResultsFolder();
  loadTwoPhaseData();
  if (drugFlowWithoutDiffusion && networkSource < 6)
    runDrugFlowWithoutDiffusion();
  if (drugFlowWithDiffusion && networkSource < 6) runDrugFlowWithDiffusion();
  if (particleFlow && networkSource < 6) runParticleFlow();
  if (angiogenesisTumour && networkSource == 4) runAngiogenesisOnLattice();
  if (angiogenesisRetina && networkSource == 5) runRetinaModel();
}

// Plotting
void network::emitPlotSignal() { emit plot(); }

// getters
pore *network::getPoreX(int i, int j, int k) const {
  if (i < 0 || i > Nx || j < 0 || j > Ny - 1 || k < 0 || k > Nz - 1) return 0;
  return tableOfPoresX[i][j][k];
}

pore *network::getPoreY(int i, int j, int k) const {
  if (i < 0 || i > Nx - 1 || j < 0 || j > Ny || k < 0 || k > Nz - 1) return 0;
  return tableOfPoresY[i][j][k];
}

pore *network::getPoreZ(int i, int j, int k) const {
  if (i < 0 || i > Nx - 1 || j < 0 || j > Ny - 1 || k < 0 || k > Nz) return 0;
  return tableOfPoresZ[i][j][k];
}

pore *network::getPoreXout(int i, int j, int k) const {
  if (i < -1 || i > Nx - 1 || j < 0 || j > Ny - 1 || k < 0 || k > Nz - 1)
    return 0;
  return tableOfPoresX[i + 1][j][k];
}

pore *network::getPoreYout(int i, int j, int k) const {
  if (i < 0 || i > Nx - 1 || j < -1 || j > Ny - 1 || k < 0 || k > Nz - 1)
    return 0;
  return tableOfPoresY[i][j + 1][k];
}

pore *network::getPoreZout(int i, int j, int k) const {
  if (i < 0 || i > Nx - 1 || j < 0 || j > Ny - 1 || k < -1 || k > Nz - 1)
    return 0;
  return tableOfPoresZ[i][j][k + 1];
}

pore *network::getPore(int i) const {
  if (i < 0 || i > totalPores - 1) return 0;
  return tableOfAllPores[i];
}

node *network::getNode(int i, int j, int k) const {
  if (i < 0 || i > Nx - 1 || j < 0 || j > Ny - 1 || k < 0 || k > Nz - 1)
    return 0;
  return tableOfNodes[i][j][k];
}

node *network::getNode(int i) const {
  if (i < 0 || i > totalNodes - 1) return 0;
  return tableOfAllNodes[i];
}

particle *network::getParticle(int i) const {
  if (i < 0 || i > totalParticles - 1) return 0;
  return tableOfParticles[i];
}

block *network::getBlock(int i, int j, int k) const {
  if (i < 0 || i > meshSizeX - 1 || j < 0 || j > meshSizeY - 1 || k < 0 ||
      k > meshSizeZ - 1)
    return 0;
  return tableOfAllBlocks[i * meshSizeY * meshSizeZ + j * meshSizeZ + k];
}

block *network::getBlock(int i) const {
  if (i < 0 || i > totalBlocks - 1) return 0;
  return tableOfAllBlocks[i];
}
int network::getTotalPores() const { return totalPores; }

int network::getTotalNodes() const { return totalNodes; }

int network::getTotalBlocks() const { return totalBlocks; }

int network::getTotalOpenedPores() const { return totalOpenedPores; }

double network::getXEdgeLength() const { return xEdgeLength; }

double network::getYEdgeLength() const { return yEdgeLength; }

double network::getZEdgeLength() const { return zEdgeLength; }

bool network::getRecord() const { return record; }

bool network::getVideoRecording() const { return videoRecording; }

int network::getNx() const { return Nx; }

void network::setNx(int value) { Nx = value; }
int network::getNy() const { return Ny; }

void network::setNy(int value) { Ny = value; }
int network::getNz() const { return Nz; }

void network::setNz(int value) { Nz = value; }

int network::getTotalParticles() const { return totalParticles; }
int network::getNetworkSource() const { return networkSource; }

bool network::getSimulationRunning() const { return simulationRunning; }

void network::setSimulationRunning(bool value) { simulationRunning = value; }
double network::getAbsolutePermeability() const { return absolutePermeability; }

double network::getPorosity() const { return porosity; }

int network::getTotalOpenedNodes() const { return totalOpenedNodes; }

bool network::getReady() const { return ready; }

void network::setCancel(bool value) { cancel = value; }
