/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "worker.h"
using namespace std;

worker::worker(network *net, int j, QObject *parent) : QObject(parent) {
  n = net;
  job = j;
}

void worker::process() {
  switch (job) {
    case 0: {
      // console mode
    }

    case 1: {
      n->setSimulationRunning(true);
      cout << "Setting up Model..." << endl;
      n->setupModel();
      cout << "Model loaded." << endl;
      n->setSimulationRunning(false);
      emit finished();
      break;
    }

    case 2: {
      n->setSimulationRunning(true);
      cout << "Starting Simulation..." << endl;
      n->runSimulation();
      cout << "End of Simulation." << endl;
      n->setSimulationRunning(false);
      emit finished();
      break;
    }
  }
}
int worker::getJob() const { return job; }

void worker::setJob(int value) { job = value; }
