#pragma once

#include <armadillo>
#include <utility>
#include <memory>

class MEl;

class IdealElementBase{
 public:
  typedef std::pair<const MEl*, const arma::mat> idealPair;
  virtual int getNumSubelements() const = 0;
  virtual const idealPair&
    getIdealPair(int subelement) const = 0;
};

class IdealElement: public IdealElementBase{
 public:
  IdealElement(const MEl* el, const arma::mat& ideal);
  int getNumSubelements() const { return 1; }
  const idealPair& getIdealPair(int subelement) const{
    return ideal_pair;
  }

 private:
  const idealPair ideal_pair;
  
};

class QuadIdealElement: public IdealElementBase{
 public:
  QuadIdealElement(const MEl* el, const arma::mat& ideal);
  int getNumSubelements() const{ return 4; }
  const idealPair& getIdealPair(int subelement) const{
    return ideal_pairs[subelement];
  }
 private:
  std::vector<idealPair> ideal_pairs;
  //idealPair ideal_pairs[4];
  std::unique_ptr<MEl> elptrs[4];
};

class PrismIdealElement: public IdealElementBase{
 public:
  PrismIdealElement(const MEl* el, const arma::mat& ideal);
  int getNumSubelements() const{ return 6; }
  const idealPair& getIdealPair(int subelement) const{
    return ideal_pairs[subelement];
  }
 private:
  std::vector<idealPair> ideal_pairs;
  //idealPair ideal_pairs[4];
  std::unique_ptr<MEl> elptrs[6];
};

