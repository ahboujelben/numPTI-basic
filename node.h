/////////////////////////////////////////////////////////////////////////////
/// Author:      Ahmed Hamdi Boujelben <ahmed.hamdi.boujelben@gmail.com>
/// Created:     2016
/// Copyright:   (c) 2020 Ahmed Hamdi Boujelben
/// Licence:     Attribution-NonCommercial 4.0 International
/////////////////////////////////////////////////////////////////////////////

#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <vector>
#include "element.h"

class node : public element {
 public:
  node(double, double, double);
  ~node();

  int getIndexX() const;
  void setIndexX(int value);

  int getIndexY() const;
  void setIndexY(int value);

  int getIndexZ() const;
  void setIndexZ(int value);

  double getXCoordinate() const;
  void setXCoordinate(double value);

  double getYCoordinate() const;
  void setYCoordinate(double value);

  double getZCoordinate() const;
  void setZCoordinate(double value);

  int getRank() const;
  void setRank(int value);

  double getPressure() const;
  void setPressure(double value);

  int getConnectionNumber() const;
  void setConnectionNumber(int value);

  std::vector<int> &getNeighboors();
  void setNeighboors(const std::vector<int> &value);

  std::vector<int> &getConnectedPores();
  void setConnectedPores(const std::vector<int> &value);

  bool getInletY() const;
  void setInletY(bool value);

  bool getOutletY() const;
  void setOutletY(bool value);

  bool getInletZ() const;
  void setInletZ(bool value);

  bool getOutletZ() const;
  void setOutletZ(bool value);

  bool getSprouTip() const;
  void setSprouTip(bool value);

  double getAge() const;
  void setAge(double value);

  node *getParent() const;
  void setParent(node *value);

  int getFeedingVesselsNumber() const;
  void setFeedingVesselsNumber(int value);

 private:
  int x;
  int y;
  int z;

  double xCoordinate;
  double yCoordinate;
  double zCoordinate;

  bool inletY;
  bool outletY;
  bool inletZ;
  bool outletZ;

  int connectionNumber;

  std::vector<int> neighboors;
  std::vector<int> connectedPores;

  int rank;
  double pressure;

  // Angiogenesis attributes
  int feedingVesselsNumber;
  bool sprouTip;
  double age;
  node *parent;
};

#endif  // NODE_H
