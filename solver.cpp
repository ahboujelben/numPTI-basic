/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"

// Eigen library
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

using namespace std;
using namespace Eigen;

void network::solvePressuresForRegularModel() {
  SparseMatrix<double> conductivityMatrix(Nx * Ny * Nz, Nx * Ny * Nz);
  conductivityMatrix.reserve(VectorXi::Constant(Nx * Ny * Nz, 7));
  VectorXd b = VectorXd::Zero(Nx * Ny * Nz);
  VectorXd pressures = VectorXd::Zero(Nx * Ny * Nz);

  int row = 0;
  for (int k = 0; k < Nz; ++k)
    for (int j = 0; j < Ny; ++j)
      for (int i = 0; i < Nx; ++i) {
        node *n, *nW, *nE, *nN, *nS, *nU, *nD;
        double c, cE, cW, cN, cS, cU, cD;
        n = getNode(i, j, k);
        nW = getNode(i - 1, j, k);
        nE = getNode(i + 1, j, k);
        nS = getNode(i, j - 1, k);
        nN = getNode(i, j + 1, k);
        nD = (Nz != 1 ? getNode(i, j, k - 1) : 0);
        nU = (Nz != 1 ? getNode(i, j, k + 1) : 0);
        cW = getPoreX(i, j, k)->getConductivity();
        cE = getPoreXout(i, j, k)->getConductivity();
        cS = getPoreY(i, j, k)->getConductivity();
        cN = getPoreYout(i, j, k)->getConductivity();
        cD = getPoreZ(i, j, k)->getConductivity();
        cU = getPoreZout(i, j, k)->getConductivity();
        c = -cW - cE - cS - cN - cD - cU;

        if (i == 0) b(row) = -pressureIn * cW;
        if (i == Nx - 1) b(row) = -pressureOut * cE;

        conductivityMatrix.insert(row, n->getRank()) = c;
        if (i != 0) conductivityMatrix.insert(row, nW->getRank()) = cW;
        if (i != Nx - 1) conductivityMatrix.insert(row, nE->getRank()) = cE;
        if (j != 0)
          if (Ny != 2) conductivityMatrix.insert(row, nS->getRank()) = cS;
        if (j != Ny - 1) conductivityMatrix.insert(row, nN->getRank()) = cN;
        if (k != 0)
          if (Nz != 1)
            if (Nz != 2) conductivityMatrix.insert(row, nD->getRank()) = cD;
        if (k != Nz - 1)
          if (Nz != 1) conductivityMatrix.insert(row, nU->getRank()) = cU;
        row++;
      }
  conductivityMatrix.makeCompressed();

  if (solverChoice == 1) {
    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(conductivityMatrix);
    pressures = solver.solve(b);
  }
  if (solverChoice == 2) {
    BiCGSTAB<SparseMatrix<double> > solver;
    solver.setTolerance(1e-6);
    solver.setMaxIterations(1000);
    solver.compute(conductivityMatrix);
    pressures = solver.solve(b);
    // cout<<solver.error()<<" "<<solver.iterations()<<endl;
  }

  for (int i = 0; i < Nx; ++i)
    for (int j = 0; j < Ny; ++j)
      for (int k = 0; k < Nz; ++k)
        getNode(i, j, k)->setPressure(pressures[getNode(i, j, k)->getRank()]);
}

void network::solvePressures() {
  SparseMatrix<double> conductivityMatrix(totalOpenedNodes, totalOpenedNodes);
  conductivityMatrix.reserve(
      VectorXi::Constant(totalOpenedNodes, maxConnectionNumber + 3));
  VectorXd b = VectorXd::Zero(totalOpenedNodes);
  VectorXd pressures = VectorXd::Zero(totalOpenedNodes);

  int row = 0;
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      vector<int> neighboors = n->getNeighboors();
      vector<int> connectedPores = n->getConnectedPores();
      double conductivity(0);
      for (unsigned j = 0; j < neighboors.size(); ++j) {
        pore* p = getPore(connectedPores[j] - 1);
        if (!p->getClosed()) {
          if (p->getInlet()) {
            b(row) += -pressureIn * p->getConductivity();
            conductivity -= p->getConductivity();
          }
          if (p->getOutlet()) {
            b(row) += -pressureOut * p->getConductivity();
            conductivity -= p->getConductivity();
          }
          if (!p->getInlet() && !p->getOutlet()) {
            node* neighboor = getNode(neighboors[j] - 1);
            if (!neighboor->getClosed()) {
              conductivityMatrix.insert(row, neighboor->getRank()) =
                  p->getConductivity();
              conductivity -= p->getConductivity();
            }
          }
        }
      }
      conductivityMatrix.insert(row, n->getRank()) = conductivity;
      row++;
    }
  }
  conductivityMatrix.makeCompressed();

  if (solverChoice == 1) {
    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(conductivityMatrix);
    pressures = solver.solve(b);
  }
  if (solverChoice == 2) {
    BiCGSTAB<SparseMatrix<double> > solver;
    solver.setTolerance(1e-6);
    solver.setMaxIterations(1000);
    solver.compute(conductivityMatrix);
    pressures = solver.solve(b);
    // cout<<solver.error()<<" "<<solver.iterations()<<endl;
  }

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) n->setPressure(pressures[n->getRank()]);
  }
}

double network::updateFlows() {
  double outletFlow(0);
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      if (p->getOutlet()) {
        p->setFlow((p->getNodeOut()->getPressure() - pressureOut) *
                   p->getConductivity());
        p->setAveragePressure((p->getNodeOut()->getPressure() + pressureOut) /
                              2);
        outletFlow += p->getFlow();
      }
      if (p->getInlet()) {
        p->setFlow((pressureIn - p->getNodeIn()->getPressure()) *
                   p->getConductivity());
        p->setAveragePressure((p->getNodeIn()->getPressure() + pressureIn) / 2);
      }
      if (!p->getInlet() && !p->getOutlet()) {
        p->setFlow(
            (p->getNodeOut()->getPressure() - p->getNodeIn()->getPressure()) *
            p->getConductivity());
        p->setAveragePressure(
            (p->getNodeIn()->getPressure() + p->getNodeOut()->getPressure()) /
            2);
      }
    }
  }
  return outletFlow;
}

void network::calculatePermeabilityAndPorosity() {
  flow = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && p->getOutlet()) flow += p->getFlow();
  }
  absolutePermeability = (flow * xEdgeLength) / (yEdgeLength * zEdgeLength *
                                                 (pressureIn - pressureOut));
  porosity = totalPoresVolume / (xEdgeLength * yEdgeLength * zEdgeLength);
  cout << "absolute permeability (mD): " << absolutePermeability / 0.987e-15
       << endl;
  cout << "porosity: " << porosity << endl;
}
