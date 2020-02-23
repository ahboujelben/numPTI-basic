/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#include <QApplication>
#include <QObject>
#include <QThread>
#include "mainwindow.h"
#include "network.h"
#include "worker.h"

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);

  MainWindow window;
  window.show();

  return a.exec();
}
