/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include "mainwindow.h"
#include <QLayout>
#include "ui_mainwindow.h"
#include "widget3d.h"

Qt::GlobalColor QtColours[] = {Qt::blue,  Qt::red,     Qt::green, Qt::gray,
                               Qt::black, Qt::magenta, Qt::yellow};

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow) {
  ui->setupUi(this);
  net = new network;
  ui->widget_3d->setNet(net);
  connect(net, SIGNAL(plot()), ui->widget_3d, SLOT(updateNetwork()));
  connect(ui->widget_3d, SIGNAL(plotted()), this, SLOT(saveImages()));
  connect(ui->widget_3d, SIGNAL(rendered()), this, SLOT(renderFinished()));
  connect(&timer, SIGNAL(timeout()), this, SLOT(plotCurvesRealTime()));
  // plot
  plotTitle = new QCPPlotTitle(ui->plotWidget, "title");
  ui->plotWidget->plotLayout()->insertRow(0);
  ui->plotWidget->plotLayout()->addElement(0, 0, plotTitle);
  ui->plotWidget->xAxis->setLabelFont(font());
  ui->plotWidget->xAxis->setTickLabelFont(font());
  ui->plotWidget->yAxis->setLabelFont(font());
  ui->plotWidget->yAxis->setTickLabelFont(font());
  ui->plotWidget->legend->setFont(font());
  ui->plotWidget->legend->setVisible(true);
  this->setWindowTitle("numPTI");

  imageIndex = 0;
  totalCurves = 0;

  timer.start(500);
}

MainWindow::~MainWindow() {
  delete ui;
  delete net;
  delete plotTitle;
}

void MainWindow::closeEvent(QCloseEvent *event) {
  if (net->getSimulationRunning()) {
    net->setCancel(true);
    while (net->getSimulationRunning()) {
      usleep(200);
    }
  }
  event->accept();
}

void MainWindow::setup3dWidget() {
  double aspect = qMax(qMax(net->getXEdgeLength(), net->getYEdgeLength()),
                       net->getZEdgeLength());
  double xInitTran = -0.5 * net->getXEdgeLength() / aspect;
  double yInitTran = -0.5 * net->getYEdgeLength() / aspect;
  double zInitTran = -0.5 * net->getZEdgeLength() / aspect;

  ui->widget_3d->setAspect(aspect);
  ui->widget_3d->setXInitTran(xInitTran);
  ui->widget_3d->setYInitTran(yInitTran);
  ui->widget_3d->setZInitTran(zInitTran);

  ui->widget_3d->setXRot(0);
  ui->widget_3d->setYRot(0);
  ui->widget_3d->setZRot(0);
  ui->widget_3d->setXTran(0);
  ui->widget_3d->setYTran(0);
  ui->widget_3d->setScale(0.4);

  if (net->getNz() == 1 && net->getNetworkSource() != 3)  //
  {
    ui->widget_3d->setXInitRot(0);
    ui->widget_3d->setYInitRot(0);
    ui->widget_3d->setScale(0.9);
  } else {
    ui->widget_3d->setXInitRot(20. / 180. * 3.14);
    ui->widget_3d->setYInitRot(150. / 180. * 3.14);
    ui->widget_3d->setXInitTran(xInitTran + 0.1);
  }
}

void MainWindow::exportNetworkDataFromGUI() {
  QSettings settings("Input Data/network_data.txt", QSettings::IniFormat);

  settings.beginGroup("Network_Source");
  int source;
  if (ui->regularRadioButton->isChecked()) source = 1;    // regular
  if (ui->fracturedRadioButton->isChecked()) source = 2;  // artificial
  if (ui->emptyNetwork->isChecked()) source = 4;   // parent vessel tumour
  if (ui->retinaNetwork->isChecked()) source = 5;  // parent vessel retina
  settings.setValue("source", source);
  settings.endGroup();

  settings.beginGroup("Geometry");
  settings.setValue("Nx", ui->NxSpinBox->text());
  settings.setValue("Ny", ui->NySpinBox->text());
  settings.setValue("Nz", ui->NzSpinBox->text());
  settings.setValue("minRadius", ui->minRadiusLineEdit->text());
  settings.setValue("maxRadius", ui->maxRadiusLineEdit->text());
  int distribution;
  if (ui->psdUniformRadioButton->isChecked()) distribution = 1;
  if (ui->psdRayleighRadioButton->isChecked()) distribution = 2;
  if (ui->psdTriangularRadioButton->isChecked()) distribution = 3;
  if (ui->psdNormalRadioButton->isChecked()) distribution = 4;
  settings.setValue("radiusDistribution", distribution);
  settings.setValue("rayleighParameter", ui->rayleighParamLineEdit->text());
  settings.setValue("triangularParameter", ui->triangularParameter->text());
  settings.setValue("normalMuParameter", ui->muNormalParameter->text());
  settings.setValue("normalSigmaParameter", ui->sigmaNormalParameter->text());
  settings.setValue("poreVolumeConstant", ui->volConsRadiusLineEdit->text());
  settings.setValue("poreVolumeExponent", ui->volExpLineEdit->text());
  settings.setValue("poreConductivityConstant", ui->condConsLineEdit->text());
  settings.setValue("poreConductivityExponent", ui->condExpLineEdit->text());
  settings.setValue("coordinationNumber", ui->zSpinBox->text());
  settings.setValue("degreeOfDistortion", ui->distortionSpinBox->text());
  settings.setValue("aspectRatio", ui->aspectRatioSpinBox->text());
  settings.setValue("length", ui->lengthLineEdit->text());
  settings.setValue("seed", ui->seed->text());
  int solverChoice;
  if (ui->choleskyRadioButton->isChecked()) solverChoice = 1;
  if (ui->bicstabRadioButton->isChecked()) solverChoice = 2;
  settings.setValue("solverChoice", solverChoice);
  settings.setValue("absolutePermeabilityCalculation",
                    ui->calcPermCheckBox->isChecked());
  settings.setValue("extractData", ui->extractDataCheckBox->isChecked());
  settings.setValue("extractionTimestep",
                    ui->extractionTimeStepSpinBox->text());
  settings.endGroup();

  settings.beginGroup("Tissue");
  settings.setValue("generateTissue", ui->generateTissue->isChecked());
  settings.setValue("meshSizeX", ui->meshSizeX->text());
  settings.setValue("meshSizeY", ui->meshSizeY->text());
  settings.setValue("meshSizeZ", ui->meshSizeZ->text());
  settings.setValue("xEdgeLength", ui->Dx->text());
  settings.setValue("yEdgeLength", ui->Dy->text());
  settings.setValue("zEdgeLength", ui->Dz->text());

  settings.endGroup();
}

void MainWindow::exportTwoPhaseDataFromGUI() {
  QSettings settings("Input Data/twoPhaseFlow_data.txt", QSettings::IniFormat);

  settings.beginGroup("Cycles");
  settings.setValue("drugFlowWithoutDiffusion",
                    ui->drugFlowWithout->isChecked());
  settings.setValue("drugFlowWithDiffusion", ui->drugFlowWith->isChecked());
  settings.setValue("particleFlow", ui->particleFlow->isChecked());
  settings.setValue("angiogenesisTumour",
                    ui->angiogenesisOnLattice->isChecked());
  settings.setValue("angiogenesisRetina",
                    ui->angiogenesisOffLattice->isChecked());
  settings.endGroup();

  settings.beginGroup("Parameters");
  settings.setValue("videoRecording",
                    ui->twoPhasePostProcessCheckBox->isChecked());
  settings.endGroup();

  settings.beginGroup("Fluids");
  settings.setValue("plasmaViscosity", ui->plasmaViscosity->text());
  settings.endGroup();

  settings.beginGroup("Flow");
  settings.setValue("flowRate", ui->twoPhaseFlowRateLineEdit->text());
  settings.setValue("simulationTime",
                    ui->twoPhaseSimulationTimeLineEdit->text());
  settings.endGroup();

  settings.beginGroup("Drug");
  settings.setValue("bolusInjection", ui->bolusRadioButton->isChecked());
  settings.setValue("bolusDuration",
                    ui->twoPhasePrecisionFactorLineEdit->text());
  settings.setValue("AIFInjection", ui->AIFRadioButton->isChecked());
  settings.endGroup();

  settings.beginGroup("Tissue");
  settings.setValue("PVT", ui->twoPhasePVT->text());
  settings.setValue("DT", ui->twoPhaseDT->text());
  settings.setValue("sigma", ui->twoPhaseSigma->text());
  settings.setValue("tissueHomo", ui->tissueHomo->isChecked());
  settings.setValue("tissueRandom", ui->tissueRandom->isChecked());
  settings.setValue("tissueCircular", ui->tissueCircular->isChecked());
  settings.setValue("tissueCircularR", ui->tissueCircularR->text());
  settings.setValue("tissueCircularX", ui->tissueCircularX->text());
  settings.setValue("tissueCircularY", ui->tissueCircularY->text());
  settings.setValue("tissueCircularZ", ui->tissueCircularZ->text());
  settings.setValue("closedBoundaries", ui->closedBoundaries->isChecked());
  settings.endGroup();

  settings.beginGroup("Particles");
  settings.setValue("injectParticles", ui->injectParticles->isChecked());
  settings.setValue("injectionInterval", ui->injectionInterval->text());
  settings.endGroup();

  settings.beginGroup("AngiogenesisOnLattice");
  settings.setValue("angio_D", ui->angio_D->text());
  settings.setValue("angio_Chi", ui->angio_Chi->text());
  settings.setValue("angio_Delta", ui->angio_Delta->text());
  settings.setValue("angio_Rho", ui->angio_Rho->text());
  settings.setValue("angio_Eta", ui->angio_Eta->text());
  settings.setValue("angio_Beta", ui->angio_Beta->text());
  settings.setValue("angio_Gamma", ui->angio_Gamma->text());
  settings.setValue("angio_Alpha", ui->angio_Alpha->text());
  settings.setValue("angio_Epsilon", ui->angio_Epsilon->text());
  settings.setValue("angio_Mu", ui->angio_Mu->text());
  settings.setValue("circularTumour", ui->circularTumour->isChecked());
  settings.setValue("linearTumour", ui->linearTumour->isChecked());
  settings.setValue("initialTipsNumber", ui->initialTipsNumber->text());
  settings.setValue("angio_Psi", ui->angio_Psi->text());
  settings.setValue("angio_TauMax", ui->angio_tauMax->text());
  settings.setValue("updateChemicals", ui->updateTissueChemicals->isChecked());
  settings.setValue("phaseSeparation", ui->phaseSeparation->isChecked());
  settings.setValue("branchingWSS", ui->branchingWSS->isChecked());
  settings.setValue("shuntPrevention", ui->shuntPrevention->isChecked());
  settings.setValue("Ks", ui->Ks->text());
  settings.setValue("Km", ui->Km->text());
  settings.setValue("Kp", ui->Kp->text());
  settings.setValue("Kc", ui->Kc->text());
  settings.setValue("Qref", ui->Qref->text());
  settings.setValue("QHDref", ui->QHDref->text());
  settings.setValue("tauRef", ui->tauRef->text());
  settings.setValue("J0", ui->J0->text());
  settings.setValue("decayConv", ui->decayConv->text());
  settings.setValue("decayCond", ui->decayCond->text());
  settings.endGroup();
}

void MainWindow::updatePlotRelativePermeabilties(string filename) {
  // prepare data
  ifstream file(filename.c_str());
  QVector<double> sw;
  QVector<double> krw;
  QVector<double> kro;
  string dummy;

  if (file) {
    getline(file, dummy);
    double value;
    while (file >> value) {
      sw.push_back(value);
      file >> value;
      krw.push_back(value);
      file >> value;
      kro.push_back(value);
    }
  }

  // remove previous graph
  ui->plotWidget->clearGraphs();

  // update title
  plotTitle->setText("Relative Permeabilities Curves");

  // create graph and assign data to it
  ui->plotWidget->addGraph();
  ui->plotWidget->graph(0)->setPen(QPen(Qt::blue));
  ui->plotWidget->graph(0)->setName("krw");
  ui->plotWidget->graph(0)->setData(sw, krw);
  ui->plotWidget->addGraph();
  ui->plotWidget->graph(1)->setPen(QPen(Qt::red));
  ui->plotWidget->graph(1)->setName("kro");
  ui->plotWidget->graph(1)->setData(sw, kro);

  // give the axes labels and ranges
  ui->plotWidget->xAxis->setLabel("Sw");
  ui->plotWidget->yAxis->setLabel("kr/k");
  ui->plotWidget->xAxis->setRange(0, 1);
  ui->plotWidget->yAxis->setRange(0, 1);

  ui->plotWidget->replot();
  ui->plotWidget->show();
}

void MainWindow::updatePlotCapillaryPressure(string filename) {
  // prepare data
  ifstream file(filename.c_str());
  QVector<double> sw;
  QVector<double> pc;
  string dummy;

  if (file) {
    getline(file, dummy);
    double value;
    while (file >> value) {
      sw.push_back(value);
      file >> value;
      pc.push_back(value);
    }
  }

  // remove previous graph
  ui->plotWidget->clearGraphs();

  // update title
  plotTitle->setText("Capillary Pressure Curve");

  // create graph and assign data to it
  ui->plotWidget->addGraph();
  ui->plotWidget->graph(0)->setData(sw, pc);
  ui->plotWidget->graph(0)->setName("pc");

  // give the axes some labels and ranges
  ui->plotWidget->xAxis->setLabel("Sw");
  ui->plotWidget->yAxis->setLabel("Pc(psi)");

  ui->plotWidget->graph(0)->rescaleAxes();
  ui->plotWidget->xAxis->setRange(0, 1);

  ui->plotWidget->replot();
  ui->plotWidget->show();
}

void MainWindow::getNetworkResults() {
  ui->permeabilityLabel->setText(
      QString::number(net->getAbsolutePermeability() / 0.987e-15, 'f', 2));
  ui->porosityLabel->setText(QString::number(net->getPorosity() * 100, 'f', 2));
  ui->totalPoresLabel->setText(QString::number(net->getTotalOpenedPores()));
  ui->totalNodesLabel->setText(QString::number(net->getTotalOpenedNodes()));
  ui->networkRunningLabel->setText("");
  ui->loadNetworkButton->setEnabled(true);
  setup3dWidget();
}

void MainWindow::getTwoPhaseSimulationResults() {
  ui->twoPhaseRunningLabel->setText("");
  ui->twoPhaseSimButton->setEnabled(true);
  ui->twoPhaseSimStopButton->setEnabled(false);
}

void MainWindow::saveImages() {
  if (net->getRecord() && net->getVideoRecording()) {
    ++imageIndex;
    QImage image = ui->widget_3d->grabFramebuffer();
    image.save("Videos/IMG" + QString::number(1000000000 + imageIndex) +
               ".png");
  }
}

void MainWindow::renderFinished() {
  ui->renderNetworkButton->setEnabled(true);
  ui->renderStopButton->setEnabled(false);
  ui->networkRenderingLabel->setText("");
  net->setSimulationRunning(false);
}

void MainWindow::on_loadNetworkButton_clicked() {
  if (net->getSimulationRunning()) return;

  ui->networkRunningLabel->setText("loading...");
  ui->loadNetworkButton->setEnabled(false);

  if (!ui->loadFromFileRadioButton->isChecked()) exportNetworkDataFromGUI();

  QThread *t = new QThread;
  worker *w = new worker(net, 1);
  w->moveToThread(t);
  connect(t, SIGNAL(started()), w, SLOT(process()));
  connect(w, SIGNAL(finished()), t, SLOT(quit()));
  connect(w, SIGNAL(finished()), this, SLOT(getNetworkResults()));
  connect(w, SIGNAL(finished()), w, SLOT(deleteLater()));
  connect(t, SIGNAL(finished()), w, SLOT(deleteLater()));
  t->start();
}

void MainWindow::on_twoPhaseSimButton_clicked() {
  if (!net->getReady() || net->getSimulationRunning()) return;

  ui->twoPhaseRunningLabel->setText("running...");
  ui->twoPhaseSimButton->setEnabled(false);
  ui->twoPhaseSimStopButton->setEnabled(true);
  imageIndex = 0;

  if (!ui->twoPhaseLoadFromFileRadioButton->isChecked())
    exportTwoPhaseDataFromGUI();

  QThread *t = new QThread;
  worker *w = new worker(net, 2);
  w->moveToThread(t);
  connect(t, SIGNAL(started()), w, SLOT(process()));
  connect(w, SIGNAL(finished()), t, SLOT(quit()));
  connect(w, SIGNAL(finished()), this, SLOT(getTwoPhaseSimulationResults()));
  connect(w, SIGNAL(finished()), w, SLOT(deleteLater()));
  connect(t, SIGNAL(finished()), w, SLOT(deleteLater()));
  t->start();
}

void MainWindow::on_twoPhaseSimStopButton_clicked() {
  net->setCancel(true);
  ui->twoPhaseRunningLabel->setText("stopping...");
}

void MainWindow::on_visibleOilCheckBox_clicked() {
  ui->widget_3d->setOilVisible(ui->visibleOilCheckBox->isChecked());
  ui->widget_3d->update();
}

void MainWindow::on_visibleWaterCheckBox_clicked() {
  ui->widget_3d->setWaterVisible(ui->visibleWaterCheckBox->isChecked());
  ui->widget_3d->update();
}

void MainWindow::on_visibleGasCheckBox_clicked() {
  ui->widget_3d->setGasVisible(ui->visibleGasCheckBox->isChecked());
  ui->widget_3d->update();
}

void MainWindow::on_pore3DCheckBox_clicked() {
  ui->widget_3d->setPoreBodies(ui->pore3DCheckBox->isChecked());
  ui->widget_3d->update();
}

void MainWindow::on_resetRadioButton_clicked() {
  setup3dWidget();
  ui->widget_3d->update();
  ui->resetRadioButton->setChecked(false);
}

void MainWindow::on_rotateCheckBox_clicked() {
  ui->widget_3d->setAnimation(ui->rotateCheckBox->isChecked());
}

void MainWindow::on_xCutCheckBox_clicked(bool checked) {
  ui->xCutSpinBox->setEnabled(checked);
  ui->xCutSpinBox_Max->setEnabled(checked);
  ui->widget_3d->setCutX(checked);
  ui->widget_3d->setCutXValue(ui->xCutSpinBox->text().toDouble());
  ui->widget_3d->setCutXValue2(ui->xCutSpinBox_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_yCutCheckBox_clicked(bool checked) {
  ui->yCutSpinBox->setEnabled(checked);
  ui->yCutSpinBox_Max->setEnabled(checked);
  ui->widget_3d->setCutY(checked);
  ui->widget_3d->setCutYValue(ui->yCutSpinBox->text().toDouble());
  ui->widget_3d->setCutYValue2(ui->yCutSpinBox_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_zCutCheckBox_clicked(bool checked) {
  ui->zCutSpinBox->setEnabled(checked);
  ui->zCutSpinBox_Max->setEnabled(checked);
  ui->widget_3d->setCutZ(checked);
  ui->widget_3d->setCutZValue(ui->zCutSpinBox->text().toDouble());
  ui->widget_3d->setCutZValue2(ui->zCutSpinBox_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_xCutSpinBox_returnPressed() {
  ui->widget_3d->setCutXValue(ui->xCutSpinBox->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_yCutSpinBox_returnPressed() {
  ui->widget_3d->setCutYValue(ui->yCutSpinBox->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_zCutSpinBox_returnPressed() {
  ui->widget_3d->setCutZValue(ui->zCutSpinBox->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_twoPhaseSaveCurveButton_clicked() {
  QString fn = QFileDialog::getSaveFileName(this, tr("Save as..."), QString(),
                                            tr("Images (*.png *.PNG)"));
  if (fn.isEmpty()) return;
  if (!(fn.endsWith(".png", Qt::CaseInsensitive))) fn += ".png";  // default

  ui->plotWidget->savePng(fn, 540, 540);
}

void MainWindow::on_saveNetworkImageButton_clicked() {
  QString fn = QFileDialog::getSaveFileName(this, tr("Save as..."), QString(),
                                            tr("Images (*.png *.PNG)"));
  if (fn.isEmpty()) return;
  if (!(fn.endsWith(".png", Qt::CaseInsensitive))) fn += ".png";  // default

  QImage image = ui->widget_3d->grabFramebuffer();
  image.save(fn, "png");
}

void MainWindow::on_renderNetworkButton_clicked() {
  if (ui->loadNetworkStatusRadioButton->isChecked()) {
    ui->widget_3d->setLoad(true);
    ui->widget_3d->setPhasePoresPath(
        ui->pathToPhasePoresLineEdit->text().toStdString());
    ui->widget_3d->setPhaseNodesPath(
        ui->pathToPhaseNodesLineEdit->text().toStdString());
    ui->widget_3d->setPhaseBlockPath(
        ui->pathToPhaseBlocksLineEdit->text().toStdString());
  }

  if (ui->renderImagesRadioButton->isChecked()) {
    ui->widget_3d->setRender(true);
    ui->widget_3d->setTotalImages(ui->imagesCountLineEdit->text().toInt());
    ui->widget_3d->setImageCount(0);
  }

  ui->renderNetworkButton->setEnabled(false);
  ui->renderStopButton->setEnabled(true);

  ui->networkRenderingLabel->setText("rendering...");
  net->setSimulationRunning(true);
}

void MainWindow::on_renderStopButton_clicked() {
  ui->widget_3d->setTotalImages(0);
  ui->networkRenderingLabel->setText("");
  net->setSimulationRunning(false);
}

void MainWindow::on_pushButton_clicked() {
  if (ui->title->text() != "") plotTitle->setText(ui->title->text());
  if (ui->xAxisTitle->text() != "")
    ui->plotWidget->xAxis->setLabel(ui->xAxisTitle->text());
  if (ui->yAxisTitle->text() != "")
    ui->plotWidget->yAxis->setLabel(ui->yAxisTitle->text());
  if (ui->minXAxis->text() != "")
    ui->plotWidget->xAxis->setRangeLower(ui->minXAxis->text().toDouble());
  if (ui->maxXAxis->text() != "")
    ui->plotWidget->xAxis->setRangeUpper(ui->maxXAxis->text().toDouble());
  if (ui->minYAxis->text() != "")
    ui->plotWidget->yAxis->setRangeLower(ui->minYAxis->text().toDouble());
  if (ui->maxYAxis->text() != "")
    ui->plotWidget->yAxis->setRangeUpper(ui->maxYAxis->text().toDouble());

  if (ui->tickStep->text() != "") {
    ui->plotWidget->xAxis->setAutoTickStep(false);
    ui->plotWidget->xAxis->setTickStep(ui->tickStep->text().toDouble());
  }
  ui->plotWidget->replot();
}

void MainWindow::on_pushButton_2_clicked() {
  ui->fileToPlot->setText(QFileDialog::getOpenFileName(
      this, tr("Select file to plot"), "Results", tr("text files (*.txt)")));
}

void MainWindow::on_pushButton_3_clicked() {
  ui->plotWidget->clearGraphs();
  ui->plotWidget->xAxis->setLabel("");
  ui->plotWidget->yAxis->setLabel("");
  plotTitle->setText("title");
  ui->plotWidget->xAxis->setAutoTickStep(true);
  ui->plotWidget->yAxis->setAutoTickStep(true);
  ui->plotWidget->rescaleAxes();
  ui->plotWidget->replot();

  totalCurves = 0;
}

void MainWindow::on_plot_clicked() {
  ifstream file(ui->fileToPlot->text().toStdString().c_str());

  if (file && totalCurves < 7) {
    string header;
    getline(file, header);

    std::stringstream ss(header);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> vstrings(begin, end);

    if (vstrings.size() <= 1) return;

    vector<QVector<double> > data;
    vector<QString> dataHeaders;
    data.resize(vstrings.size());

    for (int i = 0; i < vstrings.size(); ++i)
      dataHeaders.push_back(QString::fromStdString((vstrings[i])));

    double value;
    while (file >> value) {
      data[0].push_back(value);
      for (int i = 1; i < vstrings.size(); ++i) {
        file >> value;
        data[i].push_back(value);
      }
    }

    for (int i = 1; i < vstrings.size(); ++i) {
      QPen pen(QtColours[totalCurves]);
      pen.setWidth(2);
      ui->plotWidget->addGraph();
      ui->plotWidget->graph(totalCurves)->setPen(pen);
      ui->plotWidget->graph(totalCurves)->setName(dataHeaders[i]);
      ui->plotWidget->graph(totalCurves)->setData(data[0], data[i]);
      ui->plotWidget->rescaleAxes();
      totalCurves++;
    }
    ui->plotWidget->replot();
  }
}

void MainWindow::plotCurvesRealTime() {
  if (ui->realTimeCheckBox->isChecked()) {
    ui->plotWidget->clearGraphs();

    totalCurves = 0;

    ifstream file(ui->fileToPlot->text().toStdString().c_str());

    if (file && totalCurves < 7) {
      string header;
      getline(file, header);

      std::stringstream ss(header);
      std::istream_iterator<std::string> begin(ss);
      std::istream_iterator<std::string> end;
      std::vector<std::string> vstrings(begin, end);

      if (vstrings.size() <= 1) return;

      vector<QVector<double> > data;
      vector<QString> dataHeaders;
      data.resize(vstrings.size());

      for (int i = 0; i < vstrings.size(); ++i)
        dataHeaders.push_back(QString::fromStdString((vstrings[i])));

      double value;
      while (file >> value) {
        data[0].push_back(value);
        for (int i = 1; i < vstrings.size(); ++i) {
          file >> value;
          data[i].push_back(value);
        }
      }

      for (int i = 1; i < vstrings.size(); ++i) {
        QPen pen(QtColours[totalCurves]);
        pen.setWidth(2);
        ui->plotWidget->addGraph();
        ui->plotWidget->graph(totalCurves)->setPen(pen);
        ui->plotWidget->graph(totalCurves)->setName(dataHeaders[i]);
        ui->plotWidget->graph(totalCurves)->setData(data[0], data[i]);
        totalCurves++;
      }
      if (ui->title->text() != "") plotTitle->setText(ui->title->text());
      if (ui->xAxisTitle->text() != "")
        ui->plotWidget->xAxis->setLabel(ui->xAxisTitle->text());
      if (ui->yAxisTitle->text() != "")
        ui->plotWidget->yAxis->setLabel(ui->yAxisTitle->text());
      if (ui->minXAxis->text() != "")
        ui->plotWidget->xAxis->setRangeLower(ui->minXAxis->text().toDouble());
      if (ui->maxXAxis->text() != "")
        ui->plotWidget->xAxis->setRangeUpper(ui->maxXAxis->text().toDouble());
      if (ui->minYAxis->text() != "")
        ui->plotWidget->yAxis->setRangeLower(ui->minYAxis->text().toDouble());
      if (ui->maxYAxis->text() != "")
        ui->plotWidget->yAxis->setRangeUpper(ui->maxYAxis->text().toDouble());

      if (ui->tickStep->text() != "") {
        ui->plotWidget->xAxis->setAutoTickStep(false);
        ui->plotWidget->xAxis->setTickStep(ui->tickStep->text().toDouble());
      }
      ui->plotWidget->replot();
    }
  }
}

void MainWindow::on_xCutCheckBox_2_clicked(bool checked) {
  ui->xCutSpinBox_2->setEnabled(checked);
  ui->xCutSpinBox_2_Max->setEnabled(checked);
  ui->widget_3d->setCutXT(checked);
  ui->widget_3d->setCutXTValue(ui->xCutSpinBox_2->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_yCutCheckBox_2_clicked(bool checked) {
  ui->yCutSpinBox_2->setEnabled(checked);
  ui->yCutSpinBox_2_Max->setEnabled(checked);
  ui->widget_3d->setCutYT(checked);
  ui->widget_3d->setCutYTValue(ui->yCutSpinBox_2->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_zCutCheckBox_2_clicked(bool checked) {
  ui->zCutSpinBox_2->setEnabled(checked);
  ui->zCutSpinBox_2_Max->setEnabled(checked);
  ui->widget_3d->setCutZT(checked);
  ui->widget_3d->setCutZTValue(ui->zCutSpinBox_2->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_xCutSpinBox_2_returnPressed() {
  ui->widget_3d->setCutXTValue(ui->xCutSpinBox_2->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_yCutSpinBox_2_returnPressed() {
  ui->widget_3d->setCutYTValue(ui->yCutSpinBox_2->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_zCutSpinBox_2_returnPressed() {
  ui->widget_3d->setCutZTValue(ui->zCutSpinBox_2->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_TissueCheckBox_clicked() {
  ui->widget_3d->setTissue(ui->TissueCheckBox->isChecked());
  ui->widget_3d->update();
}

void MainWindow::on_ParticleCheckBox_clicked() {
  ui->widget_3d->setParticles(ui->ParticleCheckBox->isChecked());
  ui->widget_3d->update();
}

void MainWindow::on_xCutSpinBox_Max_returnPressed() {
  ui->widget_3d->setCutXValue2(ui->xCutSpinBox_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_yCutSpinBox_Max_returnPressed() {
  ui->widget_3d->setCutYValue2(ui->yCutSpinBox_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_zCutSpinBox_Max_returnPressed() {
  ui->widget_3d->setCutZValue2(ui->zCutSpinBox_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_xCutSpinBox_2_Max_returnPressed() {
  ui->widget_3d->setCutXTValue2(ui->xCutSpinBox_2_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_yCutSpinBox_2_Max_returnPressed() {
  ui->widget_3d->setCutYTValue2(ui->yCutSpinBox_2_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_zCutSpinBox_2_Max_returnPressed() {
  ui->widget_3d->setCutZTValue2(ui->zCutSpinBox_2_Max->text().toDouble());
  ui->widget_3d->update();
}

void MainWindow::on_visibleOilRadioButton_2_clicked() {
  ui->widget_3d->setDrugView(true);
  ui->widget_3d->setOxygenView(false);

  ui->widget_3d->setTAFView(false);
  ui->widget_3d->setFNView(false);
  ui->widget_3d->setMDEView(false);

  ui->widget_3d->update();
}

void MainWindow::on_visibleWaterRadioButton_2_clicked() {
  ui->widget_3d->setDrugView(false);
  ui->widget_3d->setOxygenView(true);

  ui->widget_3d->setTAFView(false);
  ui->widget_3d->setFNView(false);
  ui->widget_3d->setMDEView(false);

  ui->widget_3d->update();
}

void MainWindow::on_lightCheckBox_clicked(bool checked) {
  ui->widget_3d->setLight(checked);
  ui->widget_3d->update();
}

void MainWindow::on_axis3DCheckBox_clicked(bool checked) {
  ui->widget_3d->setAxis(checked);
  ui->widget_3d->update();
}

void MainWindow::on_TAF_clicked() {
  ui->widget_3d->setTAFView(true);
  ui->widget_3d->setFNView(false);
  ui->widget_3d->setMDEView(false);

  ui->widget_3d->update();
}

void MainWindow::on_FN_clicked() {
  ui->widget_3d->setTAFView(false);
  ui->widget_3d->setFNView(true);
  ui->widget_3d->setMDEView(false);

  ui->widget_3d->update();
}

void MainWindow::on_MDE_clicked() {
  ui->widget_3d->setTAFView(false);
  ui->widget_3d->setFNView(false);
  ui->widget_3d->setMDEView(true);

  ui->widget_3d->update();
}

void MainWindow::on_tumourParams_clicked() {
  ui->twoPhaseFlowRateLineEdit->setText("5e-12");
  ui->shuntPrevention->setChecked(false);
  ui->branchingWSS->setChecked(true);
  ui->Kp->setText("0.1");
  ui->Km->setText("0.07");
  ui->Kc->setText("2.6");
  ui->Ks->setText("0.35");
  ui->tauRef->setText("0.1");
  ui->Qref->setText("1e-18");
  ui->QHDref->setText("6.75e-14");
}

void MainWindow::on_retinaParams_clicked() {
  ui->twoPhaseFlowRateLineEdit->setText("15e-12");
  ui->shuntPrevention->setChecked(true);
  ui->branchingWSS->setChecked(false);
  ui->Kp->setText("0.8");
  ui->Km->setText("0.5");
  ui->Kc->setText("3.5");
  ui->Ks->setText("1");
  ui->tauRef->setText("0.5");
  ui->Qref->setText("1e-18");
  ui->QHDref->setText("15e-14");
}
