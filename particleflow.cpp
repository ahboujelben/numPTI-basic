/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

#include <thread>

using namespace std;

void network::runParticleFlow() {
  initialiseSimulation();

  ofstream file("Results/particleData.txt");
  ofstream file1("Results/vesselData.txt");
  ofstream file2("Results/nodalData.txt");
  ofstream file3("Results/exploredVolume.txt");
  file << "";

  cout << "Starting Flow in Artificial Network... " << endl;

  double startTime, endTime;
  startTime = tools::getCPUTime();

  double simulationTimeElapsed(0);
  int outputCount(0);
  double timeSoFar = 0;
  double timeToInject = 0;
  double originalStep = 100000;
  double timeToDoubleFlow = 100000;

  assignViscosities();
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

  deltaP = flowRate / A - B / A;

  pressureIn = deltaP;
  pressureOut = 0;
  solvePressures();
  updateFlows();

  // Mass conservation check

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      vector<int> connectedPores = n->getConnectedPores();
      double sumFlow = 0;
      for (unsigned j = 0; j < n->getNeighboors().size(); ++j) {
        pore* p = getPore(connectedPores[j] - 1);
        if (!p->getClosed()) {
          if (p->getNodeIn() == n)
            sumFlow += p->getFlow();
          else
            sumFlow -= p->getFlow();
        }
      }
      if (abs(sumFlow) > 1e-20) {
        cout << "Mass Not Conserved" << endl;
        cout << sumFlow << endl;
        cout << n->getId() << " " << n->getVesselType() << " " << n->getOutlet()
             << " " << n->getInlet() << endl;
        cancel = true;
        //                }
      }
    }
  }

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      vector<int> connectedPores = n->getConnectedPores();
      double sumFlowIn = 0;
      for (unsigned j = 0; j < connectedPores.size(); ++j) {
        pore* p = getPore(connectedPores[j] - 1);
        if (!p->getClosed()) {
          if ((p->getNodeIn() == n && p->getFlow() > 0) ||
              (p->getNodeOut() == n && p->getFlow() < 0)) {
            sumFlowIn += abs(p->getFlow());
          }
        }
      }
      n->setFlow(sumFlowIn);
    }
  }

  file1 << "vesselID node1ID node2ID Radius Length Flow" << endl;
  file2 << "nodeID x y z" << endl;
  file3 << "totalParticles DiscoveredVolume(%)" << endl;

  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && !p->getInlet() && !p->getOutlet()) {
      file1 << p->getId() << " " << p->getNodeIn()->getId() << " "
            << p->getNodeOut()->getId() << " " << p->getRadius() * 1e6 << " "
            << p->getLength() * 1e6 << " " << p->getFlow() << endl;
    }
  }

  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    file2 << n->getId() << " " << n->getXCoordinate() * 1e6 << " "
          << n->getYCoordinate() * 1e6 << " " << n->getZCoordinate() * 1e6
          << endl;
  }

  vector<pore*> inletPores;
  double globalTimeStep = 1e50;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && p->getInlet()) inletPores.push_back(p);
    if (!p->getClosed() && p->getVolume() / abs(p->getFlow()) < globalTimeStep)
      globalTimeStep = p->getVolume() / abs(p->getFlow());
  }

  for (int i = 0; i < totalParticles; ++i) delete tableOfParticles[i];

  totalParticles = 0;
  tableOfParticles.clear();
  double totalCurrentParticles = 0;
  double totalVisitedVolume(0);

  // post-processing
  if (videoRecording) record = true;

  while (timeSoFar < simulationTime) {
    // set time step

    timeStep = 1e50;
    for (int i = 0; i < totalParticles; ++i) {
      particle* p = getParticle(i);
      if (!p->getClosed()) {
        pore* pp = getPore(p->getPoreID() - 1);
        if (abs(pp->getFlow()) > 1e-20) {
          double fraction = pp->getFlow() > 0 ? p->getPorePosition()
                                              : 1 - p->getPorePosition();
          double step = (1 - fraction) * pp->getVolume() / abs(pp->getFlow());
          if (step < timeStep) timeStep = step;
        }
      }
    }

    timeStep = timeStep == 1e50 ? globalTimeStep : timeStep;

    if (timeStep > 1e20) cancel = true;

    // update particle positions

    for (int i = 0; i < totalParticles; ++i) {
      particle* p = getParticle(i);
      if (!p->getClosed()) {
        pore* pp = getPore(p->getPoreID() - 1);
        if (abs(pp->getFlow()) > 1e-20) {
          double newPosition =
              p->getPorePosition() + timeStep * pp->getFlow() / pp->getVolume();

          double xOut =
              pp->getNodeOut() == 0
                  ? pp->getNodeIn()->getXCoordinate() - pp->getLength()
                  : pp->getNodeOut()->getXCoordinate();
          double xIn =
              pp->getNodeIn() == 0
                  ? pp->getNodeOut()->getXCoordinate() + pp->getLength()
                  : pp->getNodeIn()->getXCoordinate();
          double yOut = pp->getNodeOut() == 0
                            ? pp->getNodeIn()->getYCoordinate()
                            : pp->getNodeOut()->getYCoordinate();
          double yIn = pp->getNodeIn() == 0 ? pp->getNodeOut()->getYCoordinate()
                                            : pp->getNodeIn()->getYCoordinate();
          double zOut = pp->getNodeOut() == 0
                            ? pp->getNodeIn()->getZCoordinate()
                            : pp->getNodeOut()->getZCoordinate();
          double zIn = pp->getNodeIn() == 0 ? pp->getNodeOut()->getZCoordinate()
                                            : pp->getNodeIn()->getZCoordinate();
          p->setXCoordinate(xOut + newPosition * (xIn - xOut));
          p->setYCoordinate(yOut + newPosition * (yIn - yOut));
          p->setZCoordinate(zOut + newPosition * (zIn - zOut));

          if (abs(newPosition) > 1e-10 && abs(newPosition) < 1 - 1e-10)
            p->setPorePosition(newPosition);

          else {
            if (!pp->getVisited()) totalVisitedVolume += pp->getVolume();
            pp->setVisited(true);
            // assign to a new pore // or close
            if (pp->getOutlet()) {
              p->setClosed(true);
              totalCurrentParticles--;
            } else {
              double dice = uniform_real();
              double cumulutativeProbability = 0;
              node* n = pp->getFlow() > 0 ? pp->getNodeIn() : pp->getNodeOut();
              vector<int> connectedPores = n->getConnectedPores();
              for (unsigned j = 0; j < connectedPores.size(); ++j) {
                pore* np = getPore(connectedPores[j] - 1);
                if (!np->getClosed()) {
                  if ((np->getNodeIn() == n && np->getFlow() < 0) ||
                      (np->getNodeOut() == n && np->getFlow() > 0)) {
                    cumulutativeProbability +=
                        abs(np->getFlow()) / n->getFlow();
                    if (dice < cumulutativeProbability) {
                      p->setPoreID(np->getId());
                      p->setPorePosition(np->getFlow() > 0 ? 0 : 1);
                      break;
                    }
                  }
                }
              }
            }
          }
        }

        if (p->getPorePosition() == 1 &&
                getPore(p->getPoreID() - 1)->getFlow() > 0 ||
            p->getPorePosition() == 0 &&
                getPore(p->getPoreID() - 1)->getFlow() < 0)
          cout << "ff" << endl;
      }
    }

    if (timeStep != 1e50) timeSoFar += timeStep;

    timeToInject += timeStep;

    if (injectParticles && timeToInject > injectionInterval) {
      totalParticles++;
      totalCurrentParticles++;
      particle* p = new particle();
      pore* pp;

      double dice = uniform_real();
      double cumulutativeProbability = 0;
      for (unsigned i = 0; i < inletPores.size(); ++i) {
        pore* inPore = inletPores[i];
        cumulutativeProbability += abs(inPore->getFlow()) / flowRate;
        if (dice < cumulutativeProbability) {
          pp = inPore;
          break;
        }
      }

      p->setPoreID(pp->getId());
      p->setPorePosition(0);
      p->setXCoordinate(pp->getNodeOut() == 0
                            ? pp->getNodeIn()->getXCoordinate() -
                                  pp->getLength()
                            : pp->getNodeOut()->getXCoordinate());
      p->setYCoordinate(pp->getNodeOut() == 0
                            ? pp->getNodeIn()->getYCoordinate()
                            : pp->getNodeOut()->getYCoordinate());
      p->setZCoordinate(pp->getNodeOut() == 0
                            ? pp->getNodeIn()->getZCoordinate()
                            : pp->getNodeOut()->getZCoordinate());
      p->setId(totalParticles);
      tableOfParticles.push_back(p);
      timeToInject = 0;

      file3 << totalParticles << " " << totalVisitedVolume / totalPoresVolume
            << endl;
    }

    emitPlotSignal();

    if (extractData && totalCurrentParticles >= 1) {
      endTime = tools::getCPUTime();
      extractParticleFlowResults(endTime - startTime, timeSoFar,
                                 simulationTimeElapsed, outputCount, true);
    }

    if (record && totalCurrentParticles >= 1)
      std::this_thread::sleep_for(std::chrono::milliseconds(5));

    if (timeSoFar > timeToDoubleFlow) {
      injectionInterval /= 2.;
      timeToDoubleFlow += originalStep;
    }

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
