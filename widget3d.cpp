/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include <libs/glew/include/GL/glew.h>
#include <boost/lexical_cast.hpp>

#include "widget3d.h"

widget3d::widget3d(QWidget *parent) : QOpenGLWidget(parent) {
  net = 0;

  light = true;
  animation = false;
  render = false;
  load = false;
  poreBodies = true;
  nodeBodies = false;
  oilVisible = true;
  waterVisible = true;
  gasVisible = true;
  waterWetVisible = true;
  oilWetVisible = true;
  cutX = false;
  cutY = false;
  cutZ = false;
  cutXT = false;
  cutYT = false;
  cutZT = false;
  particles = true;
  tissue = true;
  axis = true;

  cutXValue = 0;
  cutXValue2 = 0.5;
  cutYValue = 0;
  cutYValue2 = 0.5;
  cutZValue = 0;
  cutZValue2 = 0.5;

  cutXTValue = 0;
  cutXTValue2 = 0.5;
  cutYTValue = 0;
  cutYTValue2 = 0.5;
  cutZTValue = 0;
  cutZTValue2 = 0.5;

  drugView = true;
  oxygenView = false;
  TAFView = false;
  FNView = false;
  MDEView = false;
  connect(&timer, SIGNAL(timeout()), this, SLOT(timerUpdate()));
}

void widget3d::mousePressEvent(QMouseEvent *event) { lastPos = event->pos(); }

void widget3d::mouseMoveEvent(QMouseEvent *event) {
  int dx = event->x() - lastPos.x();
  int dy = event->y() - lastPos.y();
  if (event->buttons() & Qt::LeftButton) {
    if (QApplication::keyboardModifiers().testFlag(Qt::ControlModifier))  //
    {
      xTran += 0.002 * dx;
      yTran += -0.002 * dy;
    } else {
      xRot += 0.005 * dy;
      yRot += 0.005 * dx;
    }
    update();
  } else if (event->buttons() & Qt::RightButton) {
    xRot += 0.005 * dy;
    zRot += 0.005 * dx;
    update();
  }
  lastPos = event->pos();
}

void widget3d::wheelEvent(QWheelEvent *event) {
  double numDegrees = -event->delta() / 8.0;
  double numSteps = numDegrees / 50.0;
  scale += event->delta() / 800.;
  update();
}

void widget3d::upload(GLuint buffer, float *h_data, const unsigned count,
                      GLenum target, GLenum access) {
  glBindBuffer(target, buffer);
  glBufferData(target, count * sizeof(GLfloat), NULL, access);
  void *d_data = (void *)glMapBuffer(target, GL_READ_WRITE);
  if (d_data == NULL) {
    fprintf(stderr, "Could not map gpu buffer.\n");
    exit(1);
  }
  memcpy(d_data, (const void *)h_data, count * sizeof(GLfloat));
  if (!glUnmapBuffer(target)) {
    fprintf(stderr, "Unmap buffer failed.\n");
    exit(1);
  }
  d_data = NULL;
  glBindBuffer(target, 0);
}

unsigned widget3d::bufferVesselsData() {
  unsigned index(0);
  unsigned numberOfObjectsToDraw(0);
  if (net != 0)
    if (net->getReady()) {
      int NUMBER_CYLINDERS = net->getTotalPores();
      GLfloat *h_data = new GLfloat[11 * NUMBER_CYLINDERS];
      for (int i = 0; i < NUMBER_CYLINDERS; ++i) {
        pore *p = net->getPore(i);
        if (!p->getClosed()) {
          if (p->getNodeIn() == 0 || p->getNodeOut() == 0) continue;
          if (p->getInlet() || p->getOutlet()) continue;
          if (p->getVesselType() == 2 && !oilVisible) continue;
          if (p->getVesselType() == 1 && !waterVisible) continue;
          if (p->getVesselType() == 3 && !gasVisible) continue;
          if (cutX &&
              (p->getXCoordinate() < cutXValue * net->getXEdgeLength() ||
               p->getXCoordinate() > cutXValue2 * net->getXEdgeLength()))
            continue;
          if (cutY &&
              (p->getYCoordinate() < cutYValue * net->getYEdgeLength() ||
               p->getYCoordinate() > cutYValue2 * net->getYEdgeLength()))
            continue;
          if (cutZ &&
              (p->getZCoordinate() < cutZValue * net->getZEdgeLength() ||
               p->getZCoordinate() > cutZValue2 * net->getZEdgeLength()))
            continue;

          // center
          h_data[index] = (p->getNodeIn()->getXCoordinate() +
                           p->getNodeOut()->getXCoordinate()) /
                          2 / aspect;  // vertex.x
          h_data[index + 1] = (p->getNodeIn()->getYCoordinate() +
                               p->getNodeOut()->getYCoordinate()) /
                              2 / aspect;
          ;  // vertex.y
          h_data[index + 2] = (p->getNodeIn()->getZCoordinate() +
                               p->getNodeOut()->getZCoordinate()) /
                              2 / aspect;
          ;  // vertex.z

          // height
          h_data[index + 3] = p->getFullLength() / 2 / aspect;

          // direction
          glm::vec3 dir = glm::vec3(float((p->getNodeIn()->getXCoordinate() -
                                           p->getNodeOut()->getXCoordinate()) /
                                          aspect),
                                    float((p->getNodeIn()->getYCoordinate() -
                                           p->getNodeOut()->getYCoordinate()) /
                                          aspect),
                                    float((p->getNodeIn()->getZCoordinate() -
                                           p->getNodeOut()->getZCoordinate()) /
                                          aspect));
          h_data[index + 4] = dir[0];  // vertex.x
          h_data[index + 5] = dir[1];  // vertex.y
          h_data[index + 6] = dir[2];  // vertex.z

          // color data
          h_data[index + 7] = 1;  // phase

          float conc(0);
          if (drugView) conc = p->getConcentration();
          if (oxygenView) conc = p->getHDConcentration();
          h_data[index + 8] = conc;  // concentration1
          h_data[index + 9] = 0.0;   // concentration2
          // radius
          h_data[index + 10] = p->getRadius() / aspect;

          index += 11;
          numberOfObjectsToDraw++;
        }
      }
      if (numberOfObjectsToDraw != 0)
        upload(cylinderVBO, h_data, 11 * NUMBER_CYLINDERS, GL_ARRAY_BUFFER,
               GL_DYNAMIC_DRAW);

      delete[] h_data;
    }

  return numberOfObjectsToDraw;
}

unsigned widget3d::bufferNodesData() {
  unsigned index(0);
  unsigned numberOfObjectsToDraw(0);
  if (net != 0)
    if (net->getReady()) {
      int NUMBER_SPHERES = net->getTotalNodes();
      GLfloat *h_data = new GLfloat[7 * NUMBER_SPHERES];
      for (int i = 0; i < NUMBER_SPHERES; ++i) {
        node *p = net->getNode(i);
        if (!p->getClosed()) {
          if (p->getVesselType() == 2 && !oilVisible) continue;
          if (p->getVesselType() == 1 && !waterVisible) continue;
          if (p->getVesselType() == 3 && !gasVisible) continue;
          if (p->getRadius_sq() > pow(60e-6, 2)) continue;
          if (cutX &&
              (p->getXCoordinate() < cutXValue * net->getXEdgeLength() ||
               p->getXCoordinate() > cutXValue2 * net->getXEdgeLength()))
            continue;
          if (cutY &&
              (p->getYCoordinate() < cutYValue * net->getYEdgeLength() ||
               p->getYCoordinate() > cutYValue2 * net->getYEdgeLength()))
            continue;
          if (cutZ &&
              (p->getZCoordinate() < cutZValue * net->getZEdgeLength() ||
               p->getZCoordinate() > cutZValue2 * net->getZEdgeLength()))
            continue;

          // center
          h_data[index] = p->getXCoordinate() / aspect;      // vertex.x
          h_data[index + 1] = p->getYCoordinate() / aspect;  // vertex.y
          h_data[index + 2] = p->getZCoordinate() / aspect;  // vertex.z

          // radius
          h_data[index + 3] = p->getRadius() / aspect;

          // color data
          h_data[index + 4] = 1;  // phase

          float conc(0);
          if (drugView) conc = p->getConcentration();
          if (oxygenView) conc = p->getHDConcentration();
          h_data[index + 5] = conc;  // concentration1
          h_data[index + 6] = 0.0;   // concentration2

          index += 7;
          numberOfObjectsToDraw++;
        }
      }

      if (numberOfObjectsToDraw != 0)
        upload(sphereVBO, h_data, 7 * NUMBER_SPHERES, GL_ARRAY_BUFFER,
               GL_DYNAMIC_DRAW);

      delete[] h_data;
    }

  return numberOfObjectsToDraw;
}

unsigned widget3d::bufferBlocksData() {
  unsigned index(0);
  unsigned numberOfObjectsToDraw(0);
  if (net != 0)
    if (net->getReady()) {
      int NUMBER_CUBES = net->getTotalBlocks();
      GLfloat *h_data = new GLfloat[9 * NUMBER_CUBES];

      for (int i = 0; i < NUMBER_CUBES; ++i) {
        block *p = net->getBlock(i);
        if (!p->getClosed()) {
          if (cutXT &&
              (p->getXCoordinate() < cutXTValue * net->getXEdgeLength() ||
               p->getXCoordinate() > cutXTValue2 * net->getXEdgeLength()))
            continue;
          if (cutYT &&
              (p->getYCoordinate() < cutYTValue * net->getYEdgeLength() ||
               p->getYCoordinate() > cutYTValue2 * net->getYEdgeLength()))
            continue;
          if (cutZT &&
              (p->getZCoordinate() < cutZTValue * net->getZEdgeLength() ||
               p->getZCoordinate() > cutZTValue2 * net->getZEdgeLength()))
            continue;

          // center
          h_data[index] = p->getXCoordinate() / aspect;  // vertex.x
          h_data[index + 1] = p->getYCoordinate() / aspect;
          ;  // vertex.y
          h_data[index + 2] = p->getZCoordinate() / aspect;
          ;  // vertex.z

          // size
          h_data[index + 3] = p->getHx() / 1. / aspect;
          h_data[index + 4] = p->getHy() / 1. / aspect;
          float shrinkZ =
              (net->getNz() == 1 && net->getNetworkSource() != 3) ? 100 : 1;
          h_data[index + 5] = p->getHz() / 1. / shrinkZ / aspect;

          // color data
          h_data[index + 6] = 2;  // phase

          float conc(0);
          if (drugView) conc = p->getConcentration();

          if (oxygenView) conc = p->getHDConcentration();

          if (TAFView) conc = p->getTAFConcentration();

          if (FNView) conc = p->getFNConcentration();

          if (MDEView) conc = 1e8 * p->getMDEConcentration();

          h_data[index + 7] = conc;  // concentration1
          h_data[index + 8] = 0.0;   // concentration2

          index += 9;
          numberOfObjectsToDraw++;
        }
      }

      if (numberOfObjectsToDraw != 0)
        upload(cubeVBO, h_data, 9 * NUMBER_CUBES, GL_ARRAY_BUFFER,
               GL_DYNAMIC_DRAW);

      delete[] h_data;
    }

  return numberOfObjectsToDraw;
}

unsigned widget3d::bufferParticlesData() {
  unsigned index(0);
  unsigned numberOfObjectsToDraw(0);
  if (net != 0)
    if (net->getReady()) {
      int NUMBER_SPHERES = net->getTotalParticles();
      GLfloat *h_data = new GLfloat[7 * NUMBER_SPHERES];
      for (int i = 0; i < NUMBER_SPHERES; ++i) {
        particle *p = net->getParticle(i);
        if (!p->getClosed()) {
          // center
          h_data[index] = p->getXCoordinate() / aspect;      // vertex.x
          h_data[index + 1] = p->getYCoordinate() / aspect;  // vertex.y
          h_data[index + 2] = p->getZCoordinate() / aspect;  // vertex.z

          // radius
          h_data[index + 3] = net->getXEdgeLength() / 200. / aspect;

          // color data
          h_data[index + 4] = 3;    // phase
          h_data[index + 5] = 0.0;  // concentration1
          h_data[index + 6] = 0.0;  // concentration2

          index += 7;
          numberOfObjectsToDraw++;
        }
      }

      if (numberOfObjectsToDraw != 0)
        upload(sphereVBO, h_data, 7 * NUMBER_SPHERES, GL_ARRAY_BUFFER,
               GL_DYNAMIC_DRAW);

      delete[] h_data;
    }

  return numberOfObjectsToDraw;
}

unsigned widget3d::bufferAxesData() {
  unsigned index(0);
  unsigned numberOfObjectsToDraw(3);

  int NUMBER_CYLINDERS = 3;
  GLfloat *h_data = new GLfloat[11 * NUMBER_CYLINDERS];

  // X arrow
  // center
  h_data[index] = 0.1;    // vertex.x
  h_data[index + 1] = 0;  // vertex.y
  h_data[index + 2] = 0;  // vertex.z
  // height
  h_data[index + 3] = 0.1;
  // direction
  h_data[index + 4] = 1;  // vertex.x
  h_data[index + 5] = 0;  // vertex.y
  h_data[index + 6] = 0;  // vertex.z
  // color data
  h_data[index + 7] = 1;    // phase
  h_data[index + 8] = 0.0;  // concentration1
  h_data[index + 9] = 0.0;  // concentration2
  // radius
  h_data[index + 10] = 0.01;
  index += 11;

  // Y arrow
  // center
  h_data[index] = 0;        // vertex.x
  h_data[index + 1] = 0.1;  // vertex.y
  h_data[index + 2] = 0;    // vertex.z
  // height
  h_data[index + 3] = 0.1;
  // direction
  h_data[index + 4] = 0;  // vertex.x
  h_data[index + 5] = 1;  // vertex.y
  h_data[index + 6] = 0;  // vertex.z
  // color data
  h_data[index + 7] = 2;    // phase
  h_data[index + 8] = 0.0;  // concentration1
  h_data[index + 9] = 0.0;  // concentration2
  // radius
  h_data[index + 10] = 0.01;
  index += 11;

  // Z arrow
  // center
  h_data[index] = 0;        // vertex.x
  h_data[index + 1] = 0;    // vertex.y
  h_data[index + 2] = 0.1;  // vertex.z
  // height
  h_data[index + 3] = 0.1;
  // direction
  h_data[index + 4] = 0;  // vertex.x
  h_data[index + 5] = 0;  // vertex.y
  h_data[index + 6] = 1;  // vertex.z
  // color data
  h_data[index + 7] = 3;    // phase
  h_data[index + 8] = 0.0;  // concentration1
  h_data[index + 9] = 0.0;  // concentration2
  // radius
  h_data[index + 10] = 0.01;
  index += 11;

  upload(cylinderVBO, h_data, 11 * NUMBER_CYLINDERS, GL_ARRAY_BUFFER,
         GL_DYNAMIC_DRAW);

  delete[] h_data;

  return numberOfObjectsToDraw;
}

void widget3d::loadShaderUniforms(Shader &shader) {
  // create transformations
  glm::mat4 view(1.0f);
  glm::mat4 projection(1.0f);

  projection = glm::perspective(glm::radians(45.0f),
                                (float)width() / (float)height(), 0.1f, 100.0f);

  view = glm::translate(view, glm::vec3(0.0f, 0.0f, float(scale)));
  view = glm::translate(view, glm::vec3(0.0f + xTran, 0.0f + yTran, -2.5f));
  view = glm::rotate(view, float(xInitRot + xRot), glm::vec3(1.0f, 0.0f, 0.0f));
  view = glm::rotate(view, float(yInitRot + yRot), glm::vec3(0.0f, 1.0f, 0.0f));
  view = glm::rotate(view, 3.14f / 2.0f + float(zRot),
                     glm::vec3(0.0f, 0.0f, 1.0f));
  view = glm::translate(
      view, glm::vec3(float(xInitTran), float(yInitTran), float(zInitTran)));

  // pass transformation matrices to the shader
  shader.setMat4("projection", projection);
  shader.setMat4("view", view);

  // lighting
  shader.setVec3("lightColor", 1.0f, 1.0f, 1.0f);
  shader.setVec3("lightPos", 1.0f, 0.0f, 2.0f);

  // Normal
  shader.setVec4("eyePoint",
                 glm::inverse(view) * glm::vec4(0.0, 0.0, 0.0, 1.0));
  shader.setMat3("normalMatrix", glm::transpose(glm::inverse(glm::mat3(view))));

  // Phase Colors
  shader.setVec3("phase1Color", 0.9f, 0.35f, 0.2f);
  shader.setVec3("phase2Color", 0.4f, 0.4f, 0.5f);
  shader.setVec3("phase3Color", 0.8f, 0.8f, 0.8f);
  shader.setVec3("phase4Color", 0.65f, 0.95f, 0.15f);
}

void widget3d::loadShaderUniformsAxes(Shader &shader) {
  // create transformations
  glm::mat4 view(1.0f);
  glm::mat4 projection(1.0f);

  projection = glm::perspective(glm::radians(45.0f),
                                (float)width() / (float)height(), 0.1f, 100.0f);

  view = glm::translate(view, glm::vec3(0.0f, 0.0f, -2.5f));
  view = glm::translate(view, glm::vec3(1.0f, -1.0f, -0.5f));
  view = glm::rotate(view, float(xInitRot + xRot), glm::vec3(1.0f, 0.0f, 0.0f));
  view = glm::rotate(view, float(yInitRot + yRot), glm::vec3(0.0f, 1.0f, 0.0f));
  view = glm::rotate(view, 3.14f / 2.0f + float(zRot),
                     glm::vec3(0.0f, 0.0f, 1.0f));

  // pass transformation matrices to the shader
  shader.setMat4("projection", projection);
  shader.setMat4("view", view);

  // lighting
  shader.setVec3("lightColor", 1.0f, 1.0f, 1.0f);
  shader.setVec3("lightPos", 1.0f, 0.0f, 2.0f);

  // Normal
  shader.setVec4("eyePoint",
                 glm::inverse(view) * glm::vec4(0.0, 0.0, 0.0, 1.0));
  shader.setMat3("normalMatrix", glm::transpose(glm::inverse(glm::mat3(view))));

  // Phase Colors
  shader.setVec3("phase1Color", 0.9f, 0.2f, 0.2f);
  shader.setVec3("phase2Color", 0.2f, 0.9f, 0.4f);
  shader.setVec3("phase3Color", 0.2f, 0.2f, 0.9f);
  shader.setVec3("phase4Color", 0.7f, 1.0f, 0.15f);
}

void widget3d::drawNodes() {
  sphereShader.use();
  loadShaderUniforms(sphereShader);
  unsigned numberOfObjectsToDraw = bufferNodesData();
  glBindVertexArray(sphereVAO);
  if (numberOfObjectsToDraw != 0)
    glDrawArrays(GL_POINTS, 0, numberOfObjectsToDraw);
  glBindVertexArray(0);
}

void widget3d::drawVessels() {
  cylinderShader.use();
  loadShaderUniforms(cylinderShader);
  unsigned numberOfObjectsToDraw = bufferVesselsData();
  glBindVertexArray(cylinderVAO);
  glDrawArrays(GL_POINTS, 0, numberOfObjectsToDraw);
}

void widget3d::drawBlocks() {
  cubeShader.use();
  loadShaderUniforms(cubeShader);
  unsigned numberOfObjectsToDraw = bufferBlocksData();
  glBindVertexArray(cubeVAO);
  glDrawArrays(GL_POINTS, 0, numberOfObjectsToDraw);
}

void widget3d::drawParticles() {
  sphereShader.use();
  loadShaderUniforms(sphereShader);
  unsigned numberOfObjectsToDraw = bufferParticlesData();
  glBindVertexArray(sphereVAO);
  if (numberOfObjectsToDraw != 0)
    glDrawArrays(GL_POINTS, 0, numberOfObjectsToDraw);
  glBindVertexArray(0);
}

void widget3d::drawAxes() {
  cylinderShader.use();
  loadShaderUniformsAxes(cylinderShader);
  unsigned numberOfObjectsToDraw = bufferAxesData();
  glBindVertexArray(cylinderVAO);
  if (numberOfObjectsToDraw != 0)
    glDrawArrays(GL_POINTS, 0, numberOfObjectsToDraw);
  glBindVertexArray(0);
}

void widget3d::initializeGL() {
  glClearDepth(1.0f);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_MULTISAMPLE);

  timer.start(40);
  glewInit();
  cylinderShader.create(":/shaders/cylinder.vert", ":/shaders/cylinder.frag",
                        ":/shaders/cylinder.geom");
  sphereShader.create(":/shaders/sphere.vert", ":/shaders/sphere.frag",
                      ":/shaders/sphere.geom");
  cubeShader.create(":/shaders/cube.vert", ":/shaders/cube.frag",
                    ":/shaders/cube.geom");

  glGenVertexArrays(1, &cubeVAO);
  glGenBuffers(1, &cubeVBO);
  glBindVertexArray(cubeVAO);
  glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
  glEnableVertexAttribArray(0);  // pos
  glEnableVertexAttribArray(1);  // size
  glEnableVertexAttribArray(2);  // color
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * 4, 0);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * 4, (GLvoid *)12);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * 4, (GLvoid *)24);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  glGenVertexArrays(1, &sphereVAO);
  glGenBuffers(1, &sphereVBO);
  glBindVertexArray(sphereVAO);
  glBindBuffer(GL_ARRAY_BUFFER, sphereVBO);
  glEnableVertexAttribArray(0);  // pos
  glEnableVertexAttribArray(1);  // radius
  glEnableVertexAttribArray(2);  // color
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * 4, 0);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 7 * 4, (GLvoid *)12);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 7 * 4, (GLvoid *)16);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  glGenVertexArrays(1, &cylinderVAO);
  glGenBuffers(1, &cylinderVBO);
  glBindVertexArray(cylinderVAO);
  glBindBuffer(GL_ARRAY_BUFFER, cylinderVBO);
  glEnableVertexAttribArray(0);  // pos
  glEnableVertexAttribArray(1);  // height
  glEnableVertexAttribArray(2);  // direction
  glEnableVertexAttribArray(3);  // color
  glEnableVertexAttribArray(4);  // radius
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 11 * 4, 0);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 11 * 4, (GLvoid *)12);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 11 * 4, (GLvoid *)16);
  glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 11 * 4, (GLvoid *)28);
  glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, 11 * 4, (GLvoid *)40);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
}

void widget3d::resizeGL(int w, int h) {
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
}

void widget3d::paintGL() {
  glClearColor(16.f / 255.f, 26.f / 255.f, 32.f / 255.f, 0.0f);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  if (poreBodies) drawVessels();

  if (nodeBodies) drawNodes();

  if (tissue) drawBlocks();

  if (particles) drawParticles();

  if (axis) drawAxes();

  if (animation) yRot += 0.005;
}

bool widget3d::getLight() const { return light; }

void widget3d::setLight(bool value) { light = value; }
bool widget3d::getPoreBodies() const { return poreBodies; }

void widget3d::setPoreBodies(bool value) { poreBodies = value; }
bool widget3d::getNodeBodies() const { return nodeBodies; }

void widget3d::setNodeBodies(bool value) { nodeBodies = value; }
bool widget3d::getOilVisible() const { return oilVisible; }

void widget3d::setOilVisible(bool value) { oilVisible = value; }
bool widget3d::getWaterVisible() const { return waterVisible; }

void widget3d::setWaterVisible(bool value) { waterVisible = value; }
bool widget3d::getGasVisible() const { return gasVisible; }

void widget3d::setGasVisible(bool value) { gasVisible = value; }
bool widget3d::getWaterWetVisible() const { return waterWetVisible; }

void widget3d::setWaterWetVisible(bool value) { waterWetVisible = value; }
bool widget3d::getOilWetVisible() const { return oilWetVisible; }

void widget3d::setOilWetVisible(bool value) { oilWetVisible = value; }
double widget3d::getXRot() const { return xRot; }

void widget3d::setXRot(double value) { xRot = value; }
double widget3d::getYRot() const { return yRot; }

void widget3d::setYRot(double value) { yRot = value; }
double widget3d::getZRot() const { return zRot; }

void widget3d::setZRot(double value) { zRot = value; }
double widget3d::getXTran() const { return xTran; }

void widget3d::setXTran(double value) { xTran = value; }
double widget3d::getYTran() const { return yTran; }

void widget3d::setYTran(double value) { yTran = value; }
double widget3d::getScale() const { return scale; }

void widget3d::setScale(double value) { scale = value; }
double widget3d::getAspect() const { return aspect; }

void widget3d::setAspect(double value) { aspect = value; }

bool widget3d::getAnimation() const { return animation; }

void widget3d::setAnimation(bool value) { animation = value; }
double widget3d::getXInitTran() const { return xInitTran; }

void widget3d::setXInitTran(double value) { xInitTran = value; }
double widget3d::getYInitTran() const { return yInitTran; }

void widget3d::setYInitTran(double value) { yInitTran = value; }
double widget3d::getZInitTran() const { return zInitTran; }

void widget3d::setZInitTran(double value) { zInitTran = value; }
double widget3d::getYInitRot() const { return yInitRot; }

void widget3d::setYInitRot(double value) { yInitRot = value; }
double widget3d::getXInitRot() const { return xInitRot; }

void widget3d::setXInitRot(double value) { xInitRot = value; }
bool widget3d::getCutX() const { return cutX; }

void widget3d::setCutX(bool value) { cutX = value; }
bool widget3d::getCutY() const { return cutY; }

void widget3d::setCutY(bool value) { cutY = value; }
bool widget3d::getCutZ() const { return cutZ; }

void widget3d::setCutZ(bool value) { cutZ = value; }
double widget3d::getCutXValue() const { return cutXValue; }

void widget3d::setCutXValue(double value) { cutXValue = value; }
double widget3d::getCutYValue() const { return cutYValue; }

void widget3d::setCutYValue(double value) { cutYValue = value; }
double widget3d::getCutZValue() const { return cutZValue; }

void widget3d::setCutZValue(double value) { cutZValue = value; }

std::string widget3d::getPhasePoresPath() const { return phasePoresPath; }

void widget3d::setPhasePoresPath(const std::string &value) {
  phasePoresPath = value;
}
std::string widget3d::getPhaseNodesPath() const { return phaseNodesPath; }

void widget3d::setPhaseNodesPath(const std::string &value) {
  phaseNodesPath = value;
}
bool widget3d::getRender() const { return render; }

void widget3d::setRender(bool value) { render = value; }
bool widget3d::getLoad() const { return load; }

void widget3d::setLoad(bool value) { load = value; }
int widget3d::getTotalImages() const { return totalImages; }

void widget3d::setTotalImages(int value) { totalImages = value; }
int widget3d::getImageCount() const { return imageCount; }

void widget3d::setImageCount(int value) { imageCount = value; }

bool widget3d::getTissue() const { return tissue; }

void widget3d::setTissue(bool value) { tissue = value; }
bool widget3d::getCutXT() const { return cutXT; }

void widget3d::setCutXT(bool value) { cutXT = value; }
bool widget3d::getCutYT() const { return cutYT; }

void widget3d::setCutYT(bool value) { cutYT = value; }
bool widget3d::getCutZT() const { return cutZT; }

void widget3d::setCutZT(bool value) { cutZT = value; }
double widget3d::getCutXTValue() const { return cutXTValue; }

void widget3d::setCutXTValue(double value) { cutXTValue = value; }
double widget3d::getCutYTValue() const { return cutYTValue; }

void widget3d::setCutYTValue(double value) { cutYTValue = value; }
double widget3d::getCutZTValue() const { return cutZTValue; }

void widget3d::setCutZTValue(double value) { cutZTValue = value; }
std::string widget3d::getPhaseBlockPath() const { return phaseBlockPath; }

void widget3d::setPhaseBlockPath(const std::string &value) {
  phaseBlockPath = value;
}
bool widget3d::getParticles() const { return particles; }

void widget3d::setParticles(bool value) { particles = value; }

double widget3d::getCutXValue2() const { return cutXValue2; }

void widget3d::setCutXValue2(double value) { cutXValue2 = value; }

double widget3d::getCutYValue2() const { return cutYValue2; }

void widget3d::setCutYValue2(double value) { cutYValue2 = value; }

double widget3d::getCutZValue2() const { return cutZValue2; }

void widget3d::setCutZValue2(double value) { cutZValue2 = value; }

double widget3d::getCutXTValue2() const { return cutXTValue2; }

void widget3d::setCutXTValue2(double value) { cutXTValue2 = value; }

double widget3d::getCutYTValue2() const { return cutYTValue2; }

void widget3d::setCutYTValue2(double value) { cutYTValue2 = value; }

double widget3d::getCutZTValue2() const { return cutZTValue2; }

void widget3d::setCutZTValue2(double value) { cutZTValue2 = value; }

bool widget3d::getDrugView() const { return drugView; }

void widget3d::setDrugView(bool value) { drugView = value; }

bool widget3d::getOxygenView() const { return oxygenView; }

void widget3d::setOxygenView(bool value) { oxygenView = value; }

bool widget3d::getAxis() const { return axis; }

void widget3d::setAxis(bool value) { axis = value; }

bool widget3d::getTAFView() const { return TAFView; }

void widget3d::setTAFView(bool value) { TAFView = value; }

bool widget3d::getFNView() const { return FNView; }

void widget3d::setFNView(bool value) { FNView = value; }
bool widget3d::getMDEView() const { return MDEView; }

void widget3d::setMDEView(bool value) { MDEView = value; }

network *widget3d::getNet() const { return net; }

void widget3d::setNet(network *value) { net = value; }

void widget3d::timerUpdate() {
  if (animation || doUpdate) {
    update();
    if (net != 0)
      if (net->getReady())
        if (net->getRecord() && net->getVideoRecording()) emit plotted();
    doUpdate = false;
  }

  if (render) {
    string path = "Results/Network_Status/conc_pores" +
                  boost::lexical_cast<std::string>(1000000000 + imageCount) +
                  ".txt";
    ifstream file(path.c_str());

    if (file.is_open()) {
      int value;
      while (file >> value) {
        pore *p = net->getPore(value - 1);

        double concentration;
        file >> concentration;
        p->setConcentration(concentration);
      }

      if (net->getNetworkSource() == 2 || net->getNetworkSource() == 3) {
        string path =
            "Results/Network_Status/conc_nodes" +
            boost::lexical_cast<std::string>(1000000000 + imageCount) + ".txt";
        ifstream file(path.c_str());
        if (file.is_open()) {
          int value;
          double concentration;
          while (file >> value) {
            node *p = net->getNode(value - 1);
            file >> concentration;
            p->setConcentration(concentration);
          }
        } else {
          totalImages = 0;
        }
      }

      string path = "Results/Network_Status/conc_blocks" +
                    boost::lexical_cast<std::string>(1000000000 + imageCount) +
                    ".txt";
      ifstream file(path.c_str());
      if (file.is_open()) {
        int value;
        double conc;
        while (file >> value) {
          file >> conc;
          block *p = net->getBlock(value - 1);
          p->setConcentration(conc);
        }
      } else {
        totalImages = 0;
      }

      QImage image = grabFramebuffer();
      image.save("Videos/IMG" + QString::number(1000000000 + imageCount) +
                 ".png");

      imageCount++;
    } else {
      totalImages = 0;
    }

    if (imageCount >= totalImages) {
      render = false;
      tools::renderVideo();
      tools::cleanVideosFolder();
      emit rendered();
    }

    update();
  }

  if (load) {
    ifstream file(phasePoresPath.c_str());

    if (file.is_open()) {
      int value;
      while (file >> value) {
        pore *p = net->getPore(value - 1);
        double concentration;
        file >> concentration;
        p->setConcentration(concentration);
      }

      if (net->getNetworkSource() == 2 || net->getNetworkSource() == 3) {
        ifstream file(phaseNodesPath.c_str());
        if (file.is_open()) {
          int value;
          double concentration;
          while (file >> value) {
            node *p = net->getNode(value - 1);
            file >> concentration;
            p->setConcentration(concentration);
          }
        }
      }

      ifstream file(phaseBlockPath.c_str());
      if (file.is_open()) {
        int value;
        double conc;
        while (file >> value) {
          file >> conc;
          block *p = net->getBlock(value - 1);
          p->setConcentration(conc);
        }
      }
    }

    load = false;
    emit rendered();
    update();
  }
}

void widget3d::updateNetwork() { doUpdate = true; }
