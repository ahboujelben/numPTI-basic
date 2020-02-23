/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"
#include "tools.h"

using namespace std;

void network::runRetinaModel() {
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
