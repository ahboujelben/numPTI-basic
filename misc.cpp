/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

using namespace std;

double network::getOutletFlow() {
  double Q = 0;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && p->getOutlet()) Q += p->getFlow();
  }
  return Q;
}

// Video Extraction
void network::extractVideo() {
  if (videoRecording) {
    record = false;
    tools::renderVideo();
    tools::cleanVideosFolder();
  }
}

// randomness
int network::uniform_int(int a, int b) {
  if (a == b) return a;
  boost::random::uniform_int_distribution<> dist(a, b);
  return dist(gen);
}

double network::uniform_real(double a, double b) {
  if (a == b || a > b) return a;
  boost::random::uniform_real_distribution<> dist(a, b);
  return dist(gen);
}

double network::rayleigh(double min, double max, double ryParam) {
  if (min == max) return min;
  double value =
      min + sqrt(-pow(ryParam, 2) *
                 log(1 - uniform_real() * (1 - exp(-pow((max - min), 2) /
                                                   pow(ryParam, 2)))));
  return value;
}

double network::triangular(double a, double b, double c) {
  if (a == b || c < a || c > b) return a;
  double fc = (c - a) / (b - a);
  double u = uniform_real();
  if (u < fc)
    return a + sqrt(u * (b - a) * (c - a));
  else
    return b - sqrt((1 - u) * (b - a) * (b - c));
}

double network::normal(double min, double max, double mu, double sigma) {
  if (min == max || mu < min || mu > max) return min;

  boost::normal_distribution<> nd(mu, sigma);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
      var_nor(gen, nd);

  double value = min - 1;
  while (value < min || value > max) value = var_nor();

  return value;
}

double network::weibull(double min, double max, double alpha, double beta) {
  if (min == max) return min;

  double u = uniform_real();

  double value =
      (max - min) * pow(-beta * log(u * (1 - exp(-1 / beta)) + exp(-1 / beta)),
                        1 / alpha) +
      min;

  return value;
}

// initialisation

void network::initialiseSimulation() {
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed()) {
      p->setConcentration(0.0);
    }
  }

  for (int i = 0; i < totalNodes; ++i) {
    node* p = getNode(i);
    if (!p->getClosed()) {
      p->setConcentration(0.0);
      p->setSprouTip(false);
      p->setAge(0);
    }
  }

  for (int i = 0; i < totalBlocks; ++i) {
    block* p = getBlock(i);
    if (!p->getClosed()) {
      p->setConcentration(0.0);
      p->setHDConcentration(0.0);
      p->setMDEConcentration(0);
      p->setTAFConcentration(0);
      p->setFNConcentration(0);
    }
  }

  for (int i = 0; i < totalParticles; ++i) delete tableOfParticles[i];

  tableOfParticles.clear();
  totalParticles = 0;

  cancel = false;
}

// Data extraction
void network::extractDrugFlowResults(double timeElapsed, double timeSoFar,
                                     double& simulationTimeElapsed,
                                     int& outputCount, bool forceExtraction) {
  if (timeElapsed > simulationTimeElapsed || forceExtraction) {
    ofstream ofs1;
    ofs1.open("Results/output.txt", ofstream::app);

    ofs1 << outputCount << " " << timeSoFar << endl;

    string path = "Results/Network_Status/conc_pores" +
                  boost::lexical_cast<std::string>(1000000000 + outputCount) +
                  ".txt";

    ofstream file1(path.c_str());
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        file1 << p->getId() << " " << p->getConcentration() << endl;
      }
    }

    if (networkSource == 2 || networkSource == 3) {
      string path = "Results/Network_Status/conc_nodes" +
                    boost::lexical_cast<std::string>(1000000000 + outputCount) +
                    ".txt";
      ofstream file2(path.c_str());
      for (int i = 0; i < totalNodes; ++i) {
        node* n = getNode(i);
        if (!n->getClosed()) {
          file2 << n->getId() << " " << n->getConcentration() << endl;
        }
      }
    }

    path = "Results/Network_Status/conc_blocks" +
           boost::lexical_cast<std::string>(1000000000 + outputCount) + ".txt";

    ofstream file3(path.c_str());
    for (int i = 0; i < totalBlocks; ++i) {
      block* p = getBlock(i);
      if (!p->getClosed()) {
        file3 << p->getId() << " " << p->getConcentration() << endl;
      }
    }

    simulationTimeElapsed += extractionTimestep;
    outputCount++;

    ofs1.close();
  }
}

void network::extractParticleFlowResults(double timeElapsed, double timeSoFar,
                                         double& simulationTimeElapsed,
                                         int& outputCount,
                                         bool forceExtraction) {
  if (timeElapsed > simulationTimeElapsed || forceExtraction) {
    //        ofstream ofs1;
    //        ofs1.open("Results/output.txt", ofstream::app);

    //        ofs1<<outputCount<<" "<<timeSoFar<<endl;

    //        string
    //        path="Results/Network_Status/particles"+boost::lexical_cast<std::string>(1000000000+outputCount)+".txt";

    //        ofstream file1(path.c_str());
    //        for(int i=0;i<totalParticles;++i)
    //        {
    //            particle* p=getParticle(i);
    //            if(!p->getClosed()  )
    //            {
    //                file1<<p->getId()<<" "<<p->getXCoordinate()<<"
    //                "<<p->getYCoordinate()<<" "<<p->getZCoordinate()<<endl;
    //            }
    //        }

    //        simulationTimeElapsed+=extractionTimestep;
    //        outputCount++;

    //        ofs1.close();
    ofstream ofs1;
    ofs1.open("Results/particleData.txt", ofstream::app);

    for (int i = 0; i < totalParticles; ++i) {
      particle* p = getParticle(i);
      if (!p->getClosed()) {
        ofs1 << outputCount << "\t" << timeSoFar << "\t" << p->getId() << "\t"
             << p->getXCoordinate() * 1e6 << "\t" << p->getYCoordinate() * 1e6
             << "\t" << p->getZCoordinate() * 1e6 << endl;
      }
    }

    simulationTimeElapsed += extractionTimestep;
    outputCount++;

    ofs1.close();
  }
}

// flow rate
void network::setConstantFlowRateAker() {
  assignViscosities();

  assignConductivities();

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
}

void network::setConstantFlowRateSecant(std::set<pore*>& poresToClose) {
  assignViscosities();
  assignConductivities();

  bool stillMorePoresToClose = true;

  while (stillMorePoresToClose) {
    stillMorePoresToClose = false;

    double testQ = 0;
    double testQ1, testQ2, deltaP1, deltaP2;
    deltaP1 = deltaP * (1 + 0.5 * (flowRate - testQ) / flowRate);
    pressureIn = deltaP1;
    pressureOut = 0;
    solvePressures();
    testQ1 = updateFlows();

    deltaP2 = deltaP * (1 - 0.5 * (flowRate - testQ) / flowRate);
    pressureIn = deltaP2;
    pressureOut = 0;
    solvePressures();
    testQ2 = updateFlows();

    while (testQ > 1.01 * flowRate || testQ < 0.99 * flowRate) {
      deltaP = deltaP2 -
               (deltaP2 - deltaP1) * (testQ2 - flowRate) / (testQ2 - testQ1);
      pressureIn = deltaP;
      pressureOut = 0;
      solvePressures();
      testQ = updateFlows();
      testQ1 = testQ2;
      testQ2 = testQ;
      deltaP1 = deltaP2;
      deltaP2 = deltaP;

      // Thread Management
      if (cancel) break;
    }

    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        if (p->getOutlet() && p->getFlow() < 0) {
          p->setClosed(true);
          poresToClose.insert(p);
          stillMorePoresToClose = true;
        }
      }
    }
  }
}

void network::massConservationCheck() {
  for (int i = 0; i < totalNodes; ++i) {
    node* n = getNode(i);
    if (!n->getClosed()) {
      vector<int> connectedPores = n->getConnectedPores();
      double sumFlow = 0;
      for (unsigned j = 0; j < connectedPores.size(); ++j) {
        pore* p = getPore(connectedPores[j] - 1);
        if (!p->getClosed()) {
          if (p->getNodeIn() == n)
            sumFlow += p->getFlow();
          else
            sumFlow -= p->getFlow();
        }
      }
      if (abs(sumFlow) > 1e-13) {
        cout << "Mass Not Conserved" << endl;
        cout << sumFlow << endl;
        cancel = true;
      }
    }
  }
}
