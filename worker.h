/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef WORKER_H
#define WORKER_H

#include <QObject>
#include <iostream>
#include "network.h"

class worker : public QObject {
  Q_OBJECT
 public:
  explicit worker(network* net, int j = 0, QObject* parent = nullptr);

  int getJob() const;
  void setJob(int value);

 public slots:
  void process();

 signals:
  void finished();

 private:
  network* n;
  int job;
};

#endif  // WORKER_H
