/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef NETWORK_H
#define NETWORK_H

#include "block.h"
#include "cluster.h"
#include "node.h"
#include "particle.h"
#include "pore.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <boost/random/mersenne_twister.hpp>

#include <QObject>

using namespace std;

class network : public QObject {
  Q_OBJECT
 public:
  explicit network(QObject *parent = 0);
  ~network();
  void destroy();
  void reset();
  void setupModel();

  ////Regular Model
  void setupRegularModel();
  void createNodes();
  void createPores();
  void setNeighboors();
  void setNeighboorsForEntireNetwork();
  void applyCoordinationNumber();
  void assignRadii();
  void assignLengths();
  void distortNetwork();
  void assignShapeFactors();
  void assignShapeFactorConstants();
  void assignVolumes();
  void assignConductivities();
  void setActiveElements();
  void assignViscosities();
  void updateRanking();
  void setNeighboorsForGenericModel();
  void cleanGenericNetwork();

  ////Tissue
  void generateTissue();
  void tissueNetworkCollisionAnalysis();
  void tissueNetworkCollisionAnalysisRegular();
  void setupTissueProperties();

  ////Simulations
  void runSimulation();
  void runDrugFlowWithoutDiffusion();
  void runDrugFlowWithDiffusion();
  void runParticleFlow();
  void runAngiogenesisOnLattice();
  void runStaticModelAdaptation();
  void runRetinaModel();

  ////Mouse Network
  void setupMouseNetwork();
  void loadMouseNetworkFromFiles();

  ////Artificial Network
  void setupArtificialNetwork();
  void buildArtificialNetwork();
  void buildArtificialNetwork2();
  void buildArtificialNetwork3();
  void buildArtificialNetwork4();
  void buildArtificialNetwork5();
  void buildArtificialNetwork6();
  void buildArtificialNetwork7();

  ////Angiogenesis Network
  void createParentVessel();
  void createParentVessel3D();
  void setupTAFDistribution();
  void setupFNDistribution();
  void initialiseSroutTips(std::set<int> &);
  void calculateTimeStepForAngio();
  void updateChemicalConcentrations();
  void updateSproutTipPositions(std::set<int> &);
  void setBranching(std::set<int> &);
  void setBranchingWSS(std::set<int> &);
  void assignInitialBloodViscosities();
  void assignBloodViscosities();
  void calculateConvectedStimuli();
  void calculateConductedStimuli();
  void computeConvectiveStimuliRecursive(node *n, double stimulus);
  void computeConductiveStimuliRecursive(node *n, double stimulus);
  bool solvePressureInAngioModel();
  bool solvePressureInAngioModelWthPhaseSeparation();
  double runHaematocritFlow();
  double runHaematocritFlowWithPhaseSeparation();
  void remodelVasculature();
  bool recalculateRadii(double time = 0);
  double Chi_func(double);
  std::vector<double> generateEndothelialCellProbabilities(node *);
  node *addNode(node *, int, int, int);
  pore *addVessel(node *, node *, int);

  /// Retina Model
  void createRetinaParentVessels();

  /// Pure Diffusion/Convection
  void runTwoLayerDiffusion();
  void runLayerDissolution();

  /// Pulsating
  void pulsateVessel();

  /// Coupled Cell
  void runCoupledCell();

  ////Solvers and Permeabilities
  void solvePressures();
  void solvePressuresForRegularModel();
  double updateFlows();
  void calculatePermeabilityAndPorosity();

  ////Misc

  // tools
  double getOutletFlow();

  // Video Recording
  void extractVideo();

  // randomness
  int uniform_int(int a = 0, int b = 1);
  double uniform_real(double a = 0, double b = 1);
  double rayleigh(double, double, double);
  double triangular(double, double, double);
  double normal(double, double, double, double);
  double weibull(double, double, double, double);

  // data extraction
  void extractDrugFlowResults(double, double, double &, int &,
                              bool forceExtraction = false);
  void extractParticleFlowResults(double, double, double &, int &,
                                  bool forceExtraction = false);

  // initialisation
  void initialiseSimulation();

  // flow rate
  void setConstantFlowRateAker();
  void setConstantFlowRateSecant(std::set<pore *> &);
  void massConservationCheck();

  ////clustering
  int hkFind(int, std::vector<int> &);
  int hkUnion(std::vector<int> &, std::vector<int> &);
  int hkMakeSet(std::vector<int> &);
  // Regular
  void clusterPores(cluster *(pore::*)(void)const, void (pore::*)(cluster *),
                    char (pore::*)(void) const, char, std::vector<cluster *> &);
  void clusterEverything();

  ////loading data
  void loadData();
  void loadNetworkData();
  void loadTwoPhaseData();

  ///// Plotting
  void emitPlotSignal();

  ////Access to pores/nodes/elements
  pore *getPoreX(int, int, int) const;
  pore *getPoreY(int, int, int) const;
  pore *getPoreZ(int, int, int) const;
  pore *getPoreXout(int, int, int) const;
  pore *getPoreYout(int, int, int) const;
  pore *getPoreZout(int, int, int) const;
  pore *getPore(int) const;
  node *getNode(int, int, int) const;
  node *getNode(int) const;
  particle *getParticle(int) const;
  block *getBlock(int, int, int) const;
  block *getBlock(int) const;

  int getTotalPores() const;
  int getTotalNodes() const;
  int getTotalBlocks() const;
  int getTotalOpenedPores() const;

  ////Getters/Setters

  // ThreadManagement
  bool getReady() const;
  void setCancel(bool value);
  bool getSimulationRunning() const;
  void setSimulationRunning(bool value);

  bool getRecord() const;
  bool getVideoRecording() const;
  bool getStartRecording() const;
  void setStartRecording(bool value);

  // Getters for network attributes

  int getNetworkSource() const;
  int getTotalOpenedNodes() const;
  int getTotalParticles() const;
  double getAbsolutePermeability() const;
  double getPorosity() const;

  double getXEdgeLength() const;
  double getYEdgeLength() const;
  double getZEdgeLength() const;

  int getNx() const;
  void setNx(int value);

  int getNy() const;
  void setNy(int value);

  int getNz() const;
  void setNz(int value);

 signals:
  void plot();

 private:
  ////////////// Network Attributes //////////////
  int networkSource;
  // Basic
  int Nx;
  int Ny;
  int Nz;
  std::vector<std::vector<std::vector<node *> > > tableOfNodes;
  std::vector<std::vector<std::vector<pore *> > > tableOfPoresX;
  std::vector<std::vector<std::vector<pore *> > > tableOfPoresY;
  std::vector<std::vector<std::vector<pore *> > > tableOfPoresZ;
  std::vector<pore *> tableOfAllPores;
  std::vector<node *> tableOfAllNodes;
  std::vector<particle *> tableOfParticles;
  std::vector<block *> tableOfAllBlocks;

  int totalPores;
  int totalOpenedPores;
  int totalNodes;
  int totalOpenedNodes;
  int totalBlocks;
  int totalParticles;

  double totalPoresVolume;
  double totalNodesVolume;
  double coordinationNumber;
  double minRadius;
  double maxRadius;
  double minNodeRadius;
  double maxNodeRadius;
  int radiusDistribution;
  double length;
  double degreeOfDistortion;
  double aspectRatio;
  double shapeFactor;

  double poreVolumeConstant;
  double poreVolumeExponent;
  double poreConductivityConstant;
  double poreConductivityExponent;

  double rayleighParameter;
  double triangularParameter;
  double normalMuParameter;
  double normalSigmaParameter;

  // Extracted Network
  double xEdgeLength;
  double yEdgeLength;
  double zEdgeLength;
  int maxConnectionNumber;

  double pressureIn;
  double pressureOut;
  double flow;
  double absolutePermeability;
  double porosity;

  // seed
  int seed;

  // solver
  int solverChoice;

  // perm Calc
  bool absolutePermeabilityCalculation;

  ///////////// Flow ////////////////////////
  double flowRate;
  double deltaP;
  double timeStep;
  double simulationTime;

  ////////////// Misc Attributes //////////////
  bool record;
  bool videoRecording;
  bool extractData;
  double extractionTimestep;

  ////////////// fluids properties //////////////
  double plasmaViscosity;

  ////////////// Simulations //////////////
  bool drugFlowWithDiffusion;
  bool drugFlowWithoutDiffusion;
  bool particleFlow;
  bool angiogenesisTumour;
  bool angiogenesisRetina;
  bool cellProliferationSimulation;

  ////////////// Tissue Data //////////////
  bool buildTissue;
  double meshSizeX;
  double meshSizeY;
  double meshSizeZ;
  double PVT;
  double DT;
  double sigma;
  bool tissueHomo;
  bool tissueRandom;
  bool tissueCircular;
  double tissueCircularR;
  double tissueCircularX;
  double tissueCircularY;
  double tissueCircularZ;
  bool closedBoundaries;

  ////////////// Drug Data //////////////
  bool bolusInjection;
  double bolusDuration;
  bool AIFInjection;

  ////////////// Artificial Network //////////////
  bool injectParticles;
  double injectionInterval;

  ////////////// Angiogenesis On-Lattice data //////////////

  double h_sq;
  double angio_D;
  double angio_Chi;
  double angio_Delta;
  double angio_Rho;
  double angio_Eta;
  double angio_Beta;
  double angio_Gamma;
  double angio_Alpha;
  double angio_Epsilon;
  double angio_Mu;
  double angio_Psi;
  double angio_TauMax;
  bool circularTumour;
  bool linearTumour;
  bool updateBlockAttributes;
  bool phaseSeparation;
  bool branchingWSS;
  bool shuntPrevention;
  int initialTipsNumber;
  double Kp, Km, Ks, Kc;
  double Qref, QHDref, tauRef, J0;
  double decayConv, decayCond;

  ////////////// Clustering Attributes //////////////

  std::vector<cluster *> existClusters;

  ////////// Thread Management ///////////////
  bool cancel;
  bool ready;
  bool simulationRunning;

  ////////// Random generator ////////////////
  boost::random::mt19937 gen;
};

#endif  // NETWORK_H
