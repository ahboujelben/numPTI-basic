/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

using namespace std;

void network::runDrugFlowWithoutDiffusion() {
  initialiseSimulation();

  ofstream file("Results/output.txt");
  ofstream file1("Results/nodalPressure.txt");
  ofstream file2("Results/vesselFlows.txt");
  ofstream file3("Results/vesselConductivities.txt");
  file << "";

  ofstream ofs1, ofs2;
  ofs1.open("Results/outletConcentration.txt");
  ofs2.open("Results/averageConcentration.txt");

  ofs1 << "t OutletConc" << endl;
  ofs2 << "t AvgConc" << endl;

  cout << "Starting Flow in Network... " << endl;

  double startTime, endTime;
  startTime = tools::getCPUTime();

  double simulationTimeElapsed(0);
  int outputCount(0);
  double timeSoFar = 0;

  bool removeTracer(false);

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
      file1 << n->getId() << " " << n->getPressure() << endl;
    }
  }

  // set time step

  timeStep = 1e50;
  for (int i = 0; i < totalPores; ++i) {
    pore* p = getPore(i);
    if (!p->getClosed() && abs(p->getFlow()) > 1e-20) {
      double step = p->getVolume() / abs(p->getFlow());
      if (step < timeStep) timeStep = step;
    }

    file2 << p->getId() << " " << p->getFlow() << endl;
    file3 << p->getId() << " " << p->getConductivity() << endl;
  }

  cout << timeStep << endl;

  // post-processing
  if (videoRecording) record = true;

  int k(0);
  double oldAverageConc(0);
  double outputPV(0);
  double outputPV2(0);

  while (timeSoFar < simulationTime) {
    int fullyConcentratedPores = 0;
    int fullyConcentratedPoresBefore = 0;

    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      if (!n->getClosed()) {
        vector<int>& connectedPores = n->getConnectedPores();
        double sumFlowIn = 0;
        double massIn = 0;
        for (unsigned j = 0; j < connectedPores.size(); ++j) {
          pore* p = getPore(connectedPores[j] - 1);
          if (!p->getClosed()) {
            if ((p->getNodeIn() == n && p->getFlow() > 0) ||
                (p->getNodeOut() == n && p->getFlow() < 0)) {
              sumFlowIn += abs(p->getFlow());
              massIn += p->getConcentration() * abs(p->getFlow());
            }
          }
        }
        n->setFlow(sumFlowIn);
        n->setMassFlow(massIn);
      }
    }

    double massIn = 0;
    double flowIn = 0;

    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        if (p->getInlet())  // inlet pore
        {
          // AIF
          double A1, A2, T1, T2, sigma1, sigma2, alpha, beta, s, tau, time;
          A1 = 0.809;
          A2 = 0.33;
          T1 = 0.17046;
          T2 = 0.365;
          sigma1 = 0.0563;
          sigma2 = 0.132;
          alpha = 1.05;
          beta = 0.1685;
          s = 38.078;
          tau = 0.483;
          time = timeSoFar / 60.;

          double concentration =
              bolusInjection
                  ? 1
                  : A1 / (sigma1 * sqrt(2 * tools::pi())) *
                            exp(-pow((time - T1), 2) / (2 * pow(sigma1, 2))) +
                        A2 / (sigma2 * sqrt(2 * tools::pi())) *
                            exp(-pow((time - T2), 2) / (2 * pow(sigma2, 2))) +
                        alpha * exp(-beta * time) /
                            (1 + exp(-s * (time - tau)));
          if (AIFInjection) concentration /= 6.2;
          massIn = concentration * abs(p->getFlow());
          flowIn = abs(p->getFlow());

          if (removeTracer) flowIn = 0;
        }

        else {
          if (p->getFlow() > 0) {
            massIn = p->getNodeOut()->getMassFlow();
            flowIn = p->getNodeOut()->getFlow();
          } else {
            massIn = p->getNodeIn()->getMassFlow();
            flowIn = p->getNodeIn()->getFlow();
          }
        }

        if (abs(p->getFlow()) < 1e-20 || flowIn < 1e-20) {
          massIn = 0;
          flowIn = 1;
        }

        if (p->getConcentration() > 0.5) fullyConcentratedPoresBefore++;
        p->setConcentration(p->getConcentration() +
                            (abs(p->getFlow()) / flowIn * massIn -
                             abs(p->getFlow()) * p->getConcentration()) *
                                timeStep / p->getVolume());
        if (p->getConcentration() > 0.5) fullyConcentratedPores++;
      }
    }

    for (int i = 0; i < totalNodes; ++i) {
      node* n = getNode(i);
      if (!n->getClosed()) {
        double concentration(0);
        double neighboorsNumber(0);
        vector<int>& connectedPores = n->getConnectedPores();
        for (unsigned j = 0; j < connectedPores.size(); j++) {
          pore* npore = getPore(connectedPores[j] - 1);
          if (!npore->getClosed()) {
            concentration += npore->getConcentration();
            neighboorsNumber++;
          }
        }
        n->setConcentration(concentration / neighboorsNumber);
        // if(n->getConcentration()<-0.001 || n->getConcentration()>1.001)
        // cout<<"nfuck"<<n->getConcentration()<<endl;
      }
    }

    timeSoFar += timeStep;
    k++;

    emitPlotSignal();

    double averageConc(0), outletConc(0), totalVolume(0), totalOutletVolume(0);
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        totalVolume += p->getVolume();
        averageConc += p->getConcentration() * p->getVolume();
        if (p->getVesselType() == 3)  // vein pore
        {
          totalOutletVolume += p->getVolume();
          outletConc += p->getConcentration() * p->getVolume();
        }
      }
    }

    outletConc /= totalOutletVolume;
    averageConc /= totalVolume;

    outputPV += timeStep * flowRate;
    outputPV2 += timeStep * flowRate;
    if (outputPV > 0.04 * totalPoresVolume)
    // if((averageConc!=0 && abs(oldAverageConc-averageConc)/averageConc>0.1))
    {
      outputPV = 0;
      ofs1 << timeSoFar << " " << outletConc << endl;
      ofs2 << timeSoFar << " " << averageConc << endl;

      // oldAverageConc=averageConc;
    }

    if (extractData && outputPV2 > 0.4 * totalPoresVolume) {
      outputPV2 = 0;
      endTime = tools::getCPUTime();
      extractDrugFlowResults(endTime - startTime, timeSoFar,
                             simulationTimeElapsed, outputCount, true);
    }

    // std::this_thread::sleep_for(std::chrono::milliseconds(5));

    if (bolusInjection && timeSoFar > bolusDuration) removeTracer = true;

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

void network::runDrugFlowWithDiffusion() {
  initialiseSimulation();
  setupTissueProperties();

  ///////

  ofstream file("Results/output.txt");
  ofstream file1("Results/nodalPressure.txt");
  ofstream file2("Results/vesselFlows.txt");
  ofstream file3("Results/vesselPermeabilities.txt");
  ofstream file4("Results/blockIDs.txt");
  ofstream file5("Results/blockIDs2.txt");
  file << "";

  ofstream ofs1("Results/outletConcentration.txt");
  ofstream ofs2("Results/averageConcentration.txt");
  ofstream ofs3("Results/averageTissueConcentration.txt");
  ofstream ofs4("Results/averageVoxelConcentration.txt");

  ofs1 << "t OutletConc" << endl;
  ofs2 << "t AvgVesselConc" << endl;
  ofs3 << "t AvgTissueConc" << endl;
  ofs4 << "t AvgVoxelConc" << endl;

  cout << "Starting Flow in Network... " << endl;

  double startTime, endTime;
  startTime = tools::getCPUTime();

  double simulationTimeElapsed(0);
  int outputCount(0);
  double timeSoFar = 0;

  bool removeTracer(false);

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
      if (abs(sumFlow) > 1e-18) {
        cout << "Mass Not Conserved" << endl;
        cout << sumFlow << endl;
        cout << n->getId() << " " << n->getVesselType() << " " << n->getOutlet()
             << " " << n->getInlet() << endl;
        cancel = true;
      }
      file1 << n->getId() << " " << n->getPressure() << endl;
    }
  }

  double hx = xEdgeLength / meshSizeX;
  double hy = yEdgeLength / meshSizeY;
  double hz = zEdgeLength / meshSizeZ;

  double coefX = 1 / pow(hx, 2);
  double coefY = 1 / pow(hy, 2);
  double coefZ = 1 / pow(hz, 2);
  double coeff = 1 / pow(hx, 2) + 1 / pow(hy, 2) + 1 / pow(hz, 2);

  // set time step

  timeStep = 1e50;
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
      if (step < timeStep) timeStep = step;
    }

    file2 << p->getId() << " " << p->getFlow() << endl;
    file3 << p->getId() << " " << p->getMembranePermeability() << endl;
  }

  for (int i = 0; i < totalBlocks; ++i) {
    block* n = getBlock(i);
    file4 << n->getId() << " " << n->getXCoordinate() << " "
          << n->getYCoordinate() << " " << n->getZCoordinate() << " "
          << n->getDiffusivity() << endl;
    file5 << n->getId() << " " << n->getX() << " " << n->getY() << " "
          << n->getZ() << endl;
  }

  cout << "Time step without diffusion: " << timeStep << endl;

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
          sumSource += pp->getMembranePermeability() * area / (n->getVolume());
        }
      }

      double step = 1. / (2 * n->getDiffusivity() * coeff + sumSource + sigma);
      if (step < timeStep) timeStep = step;
    }
  }

  cout << "Time step with diffusion: " << timeStep << endl;
  double deltaT = timeStep;

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

  // post-processing
  if (videoRecording) record = true;

  int k(0);
  double outputPV(0);
  double outputPV2(0);
  while (timeSoFar < simulationTime) {
    vector<double> blockConcentration;
    vector<double> poreConcentration;
    blockConcentration.reserve(totalBlocks);
    poreConcentration.reserve(totalPores);

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
            sumSource2 += pp->getConcentration() *
                          pp->getMembranePermeability() * area /
                          (n->getVolume());
          }
        }

        double newConcentration(0);

        newConcentration += n->getConcentration() *
                            (1 - deltaT * (2 * n->getDiffusivity() * coeff +
                                           sumSource + sigma));
        if (nW != 0)
          if (!nW->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * nW->getConcentration();
        if (nE != 0)
          if (!nE->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * nE->getConcentration();
        if (nS != 0)
          if (!nS->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * nS->getConcentration();
        if (nN != 0)
          if (!nN->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * nN->getConcentration();
        if (nD != 0)
          if (!nD->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * nD->getConcentration();
        if (nU != 0)
          if (!nU->getClosed())
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * nU->getConcentration();

        if (closedBoundaries) {
          if (nW == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * n->getConcentration();
          if (nE == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefX * n->getConcentration();
          if (nS == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * n->getConcentration();
          if (nN == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefY * n->getConcentration();
          if (nD == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * n->getConcentration();
          if (nU == 0)
            newConcentration +=
                deltaT * n->getDiffusivity() * coefZ * n->getConcentration();
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
            sumSource2 += bb->getConcentration() *
                          p->getMembranePermeability() * area /
                          (p->getVolume());
          }
        }

        if (p->getInlet()) {
          if (removeTracer)
            sumInflow = 0;
          else {
            // AIF
            double A1, A2, T1, T2, sigma1, sigma2, alpha, beta, s, tau, time;
            A1 = 0.809;
            A2 = 0.33;
            T1 = 0.17046;
            T2 = 0.365;
            sigma1 = 0.0563;
            sigma2 = 0.132;
            alpha = 1.05;
            beta = 0.1685;
            s = 38.078;
            tau = 0.483;
            time = timeSoFar / 60.;

            double concentration =
                bolusInjection
                    ? 1
                    : A1 / (sigma1 * sqrt(2 * tools::pi())) *
                              exp(-pow((time - T1), 2) / (2 * pow(sigma1, 2))) +
                          A2 / (sigma2 * sqrt(2 * tools::pi())) *
                              exp(-pow((time - T2), 2) / (2 * pow(sigma2, 2))) +
                          alpha * exp(-beta * time) /
                              (1 + exp(-s * (time - tau)));
            if (AIFInjection) concentration /= 6.2;
            sumInflow = concentration * abs(p->getFlow()) / p->getVolume();
          }

        } else {
          auto& neighPores = p->getFeedingVessels();
          for (auto iterator : neighPores) {
            pore* pp = getPore(iterator.first - 1);
            double inflow = iterator.second;
            if (!pp->getClosed()) {
              sumInflow += pp->getConcentration() * p->getInflowShare() *
                           inflow / p->getVolume();
            }
          }
        }

        double newConcentration(0);
        newConcentration +=
            p->getConcentration() *
            (1 - deltaT * (abs(p->getFlow()) / p->getVolume() + sumSource));
        newConcentration += deltaT * sumInflow;
        newConcentration += deltaT * sumSource2;

        poreConcentration.push_back(newConcentration);
      }
    }

    unsigned j(0);
    for (int i = 0; i < totalBlocks; ++i) {
      block* n = getBlock(i);
      if (!n->getClosed()) {
        n->setConcentration(blockConcentration[j]);
        j++;
        if (n->getConcentration() < 0 || n->getConcentration() > 1) {
          cout << "block concentration out of range: " << n->getConcentration()
               << endl;
          cancel = true;
        }
      }
    }

    j = 0;
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        p->setConcentration(poreConcentration[j]);
        j++;
        if (p->getConcentration() < 0 || p->getConcentration() > 1.001) {
          cout << "pore concentration out of range: " << p->getConcentration()
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
            concentration += npore->getConcentration();
            neighboorsNumber++;
          }
        }
        n->setConcentration(concentration / neighboorsNumber);
        // if(n->getConcentration()<-0.001 || n->getConcentration()>1.001)
        // cout<<"nfuck"<<n->getConcentration()<<endl;
      }
    }

    timeSoFar += deltaT;
    k++;

    emitPlotSignal();

    double averageConc(0), outletConc(0), averageVoxelConc(0), totalVolume(0),
        totalOutletVolume(0);
    for (int i = 0; i < totalPores; ++i) {
      pore* p = getPore(i);
      if (!p->getClosed()) {
        totalVolume += p->getVolume();
        averageConc += p->getConcentration() * p->getVolume();
        averageVoxelConc += p->getConcentration() * p->getVolume();
        if (p->getVesselType() == 3)  // vein pore
        {
          totalOutletVolume += p->getVolume();
          outletConc += p->getConcentration() * p->getVolume();
        }
      }
    }

    double averageTissueConc(0), totalTissueVolume(0);
    for (int i = 0; i < totalBlocks; ++i) {
      block* p = getBlock(i);
      if (!p->getClosed()) {
        totalTissueVolume += p->getVolume();
        averageTissueConc += p->getConcentration() * p->getVolume();
        averageVoxelConc += p->getConcentration() * p->getVolume();
      }
    }

    outletConc /= totalOutletVolume;
    averageConc /= totalVolume;
    averageTissueConc /= totalTissueVolume;
    averageVoxelConc /= xEdgeLength * yEdgeLength * zEdgeLength;

    outputPV += timeStep * flowRate;
    outputPV2 += timeStep * flowRate;
    if (outputPV > 0.04 * totalPoresVolume)
    // if((averageConc!=0 && abs(oldAverageConc-averageConc)/averageConc>0.1) ||
    // (averageTissueConc!=0 &&
    // abs(oldTissueConc-averageTissueConc)/averageTissueConc>0.1))
    {
      outputPV = 0;
      ofs1 << timeSoFar << " " << outletConc << endl;
      ofs2 << timeSoFar << " " << averageConc << endl;
      ofs3 << timeSoFar << " " << averageTissueConc << endl;
      ofs4 << timeSoFar << " " << averageVoxelConc << endl;

      // oldAverageConc=averageConc;
      // oldTissueConc=averageTissueConc;
    }

    if (extractData && outputPV2 > 0.4 * totalPoresVolume) {
      outputPV2 = 0;
      endTime = tools::getCPUTime();
      extractDrugFlowResults(endTime - startTime, timeSoFar,
                             simulationTimeElapsed, outputCount, true);
    }

    if (bolusInjection && timeSoFar > bolusDuration) removeTracer = true;

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
