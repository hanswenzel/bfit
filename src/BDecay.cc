//#include <iostream>
#include "BDecay.h"
#include "bexclfit_fcts.h"
#include <algorithm>
#include <functional>
ClassImp(BDecay)


void BDecay::SetCuts(double c1, double c2, double c3, double c4, double c5, double c6) 
{
  min_mass = c1;
  max_mass = c2;
  chiprob = c3;
  chi2d = c4;
  pt = c5;
  ptk = c6;
}

        
void BDecay::Prune(double mmin,double mmax)
{
  Int_t flag1,flag2;
  std::cout <<"size before pruning:  "<< Size()<<std::endl;
  flag1=0; flag2=0;
  if (mmin<min_mass) {
    std::cout << "WARNING: New min. mass < current mass cut !" << std::endl;
    flag1=1;
  }
  if (mmax>max_mass) {
    std::cout << "WARNING: New max. mass > current mass cut !" << std::endl;
    flag2=1;
  }
  if (flag1 && flag2) {
    std::cout << "CANNOT APPLY CUTS ! OLD CUTS ARE MORE RESTRICT !" <<std::endl;
    return;
  }
  if (!flag1) min_mass = mmin;
  if (!flag2) max_mass = mmax;

  std::vector<Measurement>::iterator junk=remove_if(sample.begin(),sample.end(),cut(min_mass,max_mass));
  sample.erase(junk, sample.end());
  std::cout <<"size after pruning:  "<< Size()<<std::endl;

  return;
}
          
BDecay::cut::cut(double min,double max)
  : min_mass(min),max_mass(max)
{
}
      
BDecay::sbcut::sbcut(double min,double max,double masswin)
   : min_mass(min),max_mass(max),sbmasswin(masswin)
{
}     

std::vector<Measurement>* BDecay::Sidebands(double min,double max,double masswin)
{
  std::vector<Measurement> *sb = new std::vector<Measurement>;
  //sb->reserve(sample.size());
  remove_copy_if(sample.begin(),sample.end(),back_inserter(*sb),sbcut(min,max,masswin));
  std::cout << sample.size() << std::endl;
  std::cout << sb->size() << std::endl;
  return sb;
}


BDecay::~BDecay(){}

