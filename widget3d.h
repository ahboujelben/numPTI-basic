/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef WIDGET3D_H
#define WIDGET3D_H

#include "network.h"
#include "shader.h"
#include "tools.h"

#include <QApplication>
#include <QMouseEvent>
#include <QOpenGLFramebufferObject>
#include <QOpenGLWidget>
#include <QTimer>
#include <cmath>
#include <iostream>

class widget3d : public QOpenGLWidget {
  Q_OBJECT
 public:
  widget3d(QWidget *parent = 0);

 protected:
  double xRot, yRot, zRot, xTran, yTran, scale, xInitRot, yInitRot, xInitTran,
      yInitTran, zInitTran, aspect, cutXValue, cutYValue, cutZValue, cutXValue2,
      cutYValue2, cutZValue2, cutXTValue, cutYTValue, cutZTValue, cutXTValue2,
      cutYTValue2, cutZTValue2;
  bool light, animation, doUpdate, render, load, poreBodies, nodeBodies, tissue,
      particles, axis, oilVisible, waterVisible, gasVisible, waterWetVisible,
      oilWetVisible, cutX, cutY, cutZ, cutXT, cutYT, cutZT, drugView,
      oxygenView, TAFView, FNView, MDEView;
  int totalImages, imageCount;
  std::string phasePoresPath, phaseNodesPath, phaseBlockPath;
  QPoint lastPos;
  network *net;
  QTimer timer;
  QOpenGLFramebufferObject *fbo;
  unsigned int cubeVBO, cubeVAO, sphereVBO, sphereVAO, cylinderVBO, cylinderVAO;
  Shader cubeShader, sphereShader, cylinderShader;
  std::vector<float> cubeData, sphereData, cylinderData;

  void mousePressEvent(QMouseEvent *event);

  void mouseMoveEvent(QMouseEvent *event);

  void wheelEvent(QWheelEvent *event);

  void upload(GLuint buffer, float *h_data, const unsigned count, GLenum target,
              GLenum access);

  unsigned bufferVesselsData();

  unsigned bufferNodesData();

  unsigned bufferBlocksData();

  unsigned bufferParticlesData();

  unsigned bufferAxesData();

  void loadShaderUniforms(Shader &shader);

  void loadShaderUniformsAxes(Shader &shader);

  void drawNodes();

  void drawVessels();

  void drawBlocks();

  void drawParticles();

  void drawAxes();

  void initializeGL();

  void resizeGL(int w, int h);

  void paintGL();

 public:
  bool getLight() const;
  void setLight(bool value);

  bool getPoreBodies() const;
  void setPoreBodies(bool value);

  bool getNodeBodies() const;
  void setNodeBodies(bool value);

  bool getOilVisible() const;
  void setOilVisible(bool value);

  bool getWaterVisible() const;
  void setWaterVisible(bool value);

  bool getGasVisible() const;
  void setGasVisible(bool value);

  bool getWaterWetVisible() const;
  void setWaterWetVisible(bool value);

  bool getOilWetVisible() const;
  void setOilWetVisible(bool value);

  double getXRot() const;
  void setXRot(double value);

  double getYRot() const;
  void setYRot(double value);

  double getZRot() const;
  void setZRot(double value);

  double getXTran() const;
  void setXTran(double value);

  double getYTran() const;
  void setYTran(double value);

  double getScale() const;
  void setScale(double value);

  double getAspect() const;
  void setAspect(double value);

  bool getAnimation() const;
  void setAnimation(bool value);

  double getXInitTran() const;
  void setXInitTran(double value);

  double getYInitTran() const;
  void setYInitTran(double value);

  double getZInitTran() const;
  void setZInitTran(double value);

  double getYInitRot() const;
  void setYInitRot(double value);

  double getXInitRot() const;
  void setXInitRot(double value);

  bool getCutX() const;
  void setCutX(bool value);

  bool getCutY() const;
  void setCutY(bool value);

  bool getCutZ() const;
  void setCutZ(bool value);

  double getCutXValue() const;
  void setCutXValue(double value);

  double getCutYValue() const;
  void setCutYValue(double value);

  double getCutZValue() const;
  void setCutZValue(double value);

  std::string getPhasePoresPath() const;
  void setPhasePoresPath(const std::string &value);

  std::string getPhaseNodesPath() const;
  void setPhaseNodesPath(const std::string &value);

  bool getRender() const;
  void setRender(bool value);

  bool getLoad() const;
  void setLoad(bool value);

  int getTotalImages() const;
  void setTotalImages(int value);

  int getImageCount() const;
  void setImageCount(int value);

  bool getTissue() const;
  void setTissue(bool value);

  bool getCutXT() const;
  void setCutXT(bool value);

  bool getCutYT() const;
  void setCutYT(bool value);

  bool getCutZT() const;
  void setCutZT(bool value);

  double getCutXTValue() const;
  void setCutXTValue(double value);

  double getCutYTValue() const;
  void setCutYTValue(double value);

  double getCutZTValue() const;
  void setCutZTValue(double value);

  std::string getPhaseBlockPath() const;
  void setPhaseBlockPath(const std::string &value);

  bool getParticles() const;
  void setParticles(bool value);

  double getCutXValue2() const;
  void setCutXValue2(double value);

  double getCutYValue2() const;
  void setCutYValue2(double value);

  double getCutZValue2() const;
  void setCutZValue2(double value);

  double getCutXTValue2() const;
  void setCutXTValue2(double value);

  double getCutYTValue2() const;
  void setCutYTValue2(double value);

  double getCutZTValue2() const;
  void setCutZTValue2(double value);

  bool getDrugView() const;
  void setDrugView(bool value);

  bool getOxygenView() const;
  void setOxygenView(bool value);

  bool getAxis() const;
  void setAxis(bool value);

  bool getTAFView() const;
  void setTAFView(bool value);

  bool getFNView() const;
  void setFNView(bool value);

  bool getMDEView() const;
  void setMDEView(bool value);

  network *getNet() const;
  void setNet(network *value);

 public slots:
  void timerUpdate();

  void updateNetwork();

 signals:
  void plotted();
  void rendered();
};

#endif  // WIDGET3D_H
