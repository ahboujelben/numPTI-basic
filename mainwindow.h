/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "libs/qcustomplot/qcustomplot.h"
#include "network.h"
#include "worker.h"

#include <unistd.h>
#include <QMainWindow>
#include "QThread"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
  Q_OBJECT

 public:
  explicit MainWindow(QWidget *parent = nullptr);
  ~MainWindow();

  void closeEvent(QCloseEvent *event);
  void setup3dWidget();
  void exportNetworkDataFromGUI();
  void exportTwoPhaseDataFromGUI();
  void exportDepletionDataFromGUI();
  void updatePlotRelativePermeabilties(std::string);
  void updatePlotCapillaryPressure(std::string);

 private:
  Ui::MainWindow *ui;
  network *net;
  QCPPlotTitle *plotTitle;
  QTimer timer;
  int imageIndex;
  int totalCurves;

 signals:
  void closing();

 private slots:
  void getNetworkResults();
  void getTwoPhaseSimulationResults();
  void saveImages();
  void renderFinished();
  void on_loadNetworkButton_clicked();
  void on_twoPhaseSimButton_clicked();
  void on_twoPhaseSimStopButton_clicked();
  void on_visibleOilCheckBox_clicked();
  void on_visibleWaterCheckBox_clicked();
  void on_visibleGasCheckBox_clicked();
  void on_pore3DCheckBox_clicked();
  void on_resetRadioButton_clicked();
  void on_rotateCheckBox_clicked();
  void on_xCutCheckBox_clicked(bool checked);
  void on_yCutCheckBox_clicked(bool checked);
  void on_zCutCheckBox_clicked(bool checked);
  void on_xCutSpinBox_returnPressed();
  void on_yCutSpinBox_returnPressed();
  void on_zCutSpinBox_returnPressed();

  void on_twoPhaseSaveCurveButton_clicked();
  void on_saveNetworkImageButton_clicked();
  void on_renderNetworkButton_clicked();
  void on_renderStopButton_clicked();
  void on_pushButton_clicked();
  void on_pushButton_2_clicked();
  void on_pushButton_3_clicked();
  void on_plot_clicked();
  void plotCurvesRealTime();
  void on_xCutCheckBox_2_clicked(bool checked);
  void on_yCutCheckBox_2_clicked(bool checked);
  void on_zCutCheckBox_2_clicked(bool checked);
  void on_xCutSpinBox_2_returnPressed();
  void on_yCutSpinBox_2_returnPressed();
  void on_zCutSpinBox_2_returnPressed();
  void on_TissueCheckBox_clicked();
  void on_ParticleCheckBox_clicked();
  void on_xCutSpinBox_Max_returnPressed();
  void on_yCutSpinBox_Max_returnPressed();
  void on_zCutSpinBox_Max_returnPressed();
  void on_xCutSpinBox_2_Max_returnPressed();
  void on_yCutSpinBox_2_Max_returnPressed();
  void on_zCutSpinBox_2_Max_returnPressed();
  void on_visibleOilRadioButton_2_clicked();
  void on_visibleWaterRadioButton_2_clicked();
  void on_lightCheckBox_clicked(bool checked);
  void on_axis3DCheckBox_clicked(bool checked);
  void on_TAF_clicked();
  void on_FN_clicked();
  void on_MDE_clicked();
  void on_tumourParams_clicked();
  void on_retinaParams_clicked();
};

#endif  // MAINWINDOW_H
