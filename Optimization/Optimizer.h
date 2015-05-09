#pragma once

class MeritEvaluator;

class Optimizer{
 private:
 
 protected:
  MeritEvaluator& merit_evaluator;
  double merit;

 public:
 Optimizer(MeritEvaluator& merit_evaluator_t): 
  merit_evaluator(merit_evaluator_t){}
  virtual int Optimize() = 0;
  virtual double getMerit(){ return merit; }
  
};

