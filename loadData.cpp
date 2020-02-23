/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "network.h"

#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace std;

void network::loadNetworkData() {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("Input Data/network_data.txt", pt);

  networkSource = pt.get<int>("Network_Source.source");

  Nx = pt.get<int>("Geometry.Nx");
  Ny = pt.get<int>("Geometry.Ny");
  Nz = pt.get<int>("Geometry.Nz");
  minRadius = pt.get<double>("Geometry.minRadius") * 1e-6;
  maxRadius = pt.get<double>("Geometry.maxRadius") * 1e-6;
  radiusDistribution = pt.get<int>("Geometry.radiusDistribution");
  rayleighParameter = pt.get<double>("Geometry.rayleighParameter") * 1e-6;
  triangularParameter = pt.get<double>("Geometry.triangularParameter") * 1e-6;
  normalMuParameter = pt.get<double>("Geometry.normalMuParameter");
  normalSigmaParameter = pt.get<double>("Geometry.normalSigmaParameter");
  poreVolumeConstant = pt.get<double>("Geometry.poreVolumeConstant");
  poreVolumeExponent = pt.get<double>("Geometry.poreVolumeExponent");
  poreConductivityConstant =
      pt.get<double>("Geometry.poreConductivityConstant");
  poreConductivityExponent =
      pt.get<double>("Geometry.poreConductivityExponent");
  coordinationNumber = pt.get<double>("Geometry.coordinationNumber");
  degreeOfDistortion = pt.get<double>("Geometry.degreeOfDistortion");
  aspectRatio = pt.get<double>("Geometry.aspectRatio");
  length = pt.get<double>("Geometry.length") * 1e-6;
  seed = pt.get<int>("Geometry.seed");
  solverChoice = pt.get<int>("Geometry.solverChoice");
  absolutePermeabilityCalculation =
      pt.get<bool>("Geometry.absolutePermeabilityCalculation");
  extractData = pt.get<bool>("Geometry.extractData");
  extractionTimestep = pt.get<double>("Geometry.extractionTimestep");

  buildTissue = pt.get<bool>("Tissue.generateTissue");
  meshSizeX = pt.get<double>("Tissue.meshSizeX");
  meshSizeY = pt.get<double>("Tissue.meshSizeY");
  meshSizeZ = pt.get<double>("Tissue.meshSizeZ");
  xEdgeLength = pt.get<double>("Tissue.xEdgeLength") * 1e-6;
  yEdgeLength = pt.get<double>("Tissue.yEdgeLength") * 1e-6;
  zEdgeLength = pt.get<double>("Tissue.zEdgeLength") * 1e-6;
}

void network::loadTwoPhaseData() {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini("Input Data/twoPhaseFlow_data.txt",
                                             pt);

  drugFlowWithoutDiffusion = pt.get<bool>("Cycles.drugFlowWithoutDiffusion");
  drugFlowWithDiffusion = pt.get<bool>("Cycles.drugFlowWithDiffusion");
  particleFlow = pt.get<bool>("Cycles.particleFlow");
  angiogenesisTumour = pt.get<bool>("Cycles.angiogenesisTumour");
  angiogenesisRetina = pt.get<bool>("Cycles.angiogenesisRetina");

  videoRecording = pt.get<bool>("Parameters.videoRecording");

  plasmaViscosity = pt.get<double>("Fluids.plasmaViscosity") * 1e-3;

  flowRate = pt.get<double>("Flow.flowRate");
  simulationTime = pt.get<double>("Flow.simulationTime");

  bolusDuration = pt.get<double>("Drug.bolusDuration");
  bolusInjection = pt.get<bool>("Drug.bolusInjection");
  AIFInjection = pt.get<bool>("Drug.AIFInjection");

  PVT = pt.get<double>("Tissue.PVT");
  DT = pt.get<double>("Tissue.DT");
  sigma = pt.get<double>("Tissue.sigma");
  tissueHomo = pt.get<bool>("Tissue.tissueHomo");
  tissueRandom = pt.get<bool>("Tissue.tissueRandom");
  tissueCircular = pt.get<bool>("Tissue.tissueCircular");
  tissueCircularR = pt.get<double>("Tissue.tissueCircularR") * 1e-6;
  tissueCircularX = pt.get<double>("Tissue.tissueCircularX");
  tissueCircularY = pt.get<double>("Tissue.tissueCircularY");
  tissueCircularZ = pt.get<double>("Tissue.tissueCircularZ");
  closedBoundaries = pt.get<bool>("Tissue.closedBoundaries");

  injectParticles = pt.get<bool>("Particles.injectParticles");
  injectionInterval = pt.get<double>("Particles.injectionInterval");

  angio_D = pt.get<double>("AngiogenesisOnLattice.angio_D");
  angio_Chi = pt.get<double>("AngiogenesisOnLattice.angio_Chi");
  angio_Delta = pt.get<double>("AngiogenesisOnLattice.angio_Delta");
  angio_Rho = pt.get<double>("AngiogenesisOnLattice.angio_Rho");
  angio_Eta = pt.get<double>("AngiogenesisOnLattice.angio_Eta");
  angio_Beta = pt.get<double>("AngiogenesisOnLattice.angio_Beta");
  angio_Gamma = pt.get<double>("AngiogenesisOnLattice.angio_Gamma");
  angio_Alpha = pt.get<double>("AngiogenesisOnLattice.angio_Alpha");
  angio_Epsilon = pt.get<double>("AngiogenesisOnLattice.angio_Epsilon");
  angio_Mu = pt.get<double>("AngiogenesisOnLattice.angio_Mu");
  angio_Psi = pt.get<double>("AngiogenesisOnLattice.angio_Psi");
  angio_TauMax = pt.get<double>("AngiogenesisOnLattice.angio_TauMax");
  circularTumour = pt.get<bool>("AngiogenesisOnLattice.circularTumour");
  linearTumour = pt.get<bool>("AngiogenesisOnLattice.linearTumour");
  updateBlockAttributes = pt.get<bool>("AngiogenesisOnLattice.updateChemicals");
  phaseSeparation = pt.get<bool>("AngiogenesisOnLattice.phaseSeparation");
  branchingWSS = pt.get<bool>("AngiogenesisOnLattice.branchingWSS");
  shuntPrevention = pt.get<bool>("AngiogenesisOnLattice.shuntPrevention");
  initialTipsNumber = pt.get<int>("AngiogenesisOnLattice.initialTipsNumber");
  Kp = pt.get<double>("AngiogenesisOnLattice.Kp");
  Ks = pt.get<double>("AngiogenesisOnLattice.Ks");
  Km = pt.get<double>("AngiogenesisOnLattice.Km");
  Kc = pt.get<double>("AngiogenesisOnLattice.Kc");
  Qref = pt.get<double>("AngiogenesisOnLattice.Qref");
  QHDref = pt.get<double>("AngiogenesisOnLattice.QHDref");
  tauRef = pt.get<double>("AngiogenesisOnLattice.tauRef");
  J0 = pt.get<double>("AngiogenesisOnLattice.J0");
  decayConv = pt.get<double>("AngiogenesisOnLattice.decayConv");
  decayCond = pt.get<double>("AngiogenesisOnLattice.decayCond");
}
