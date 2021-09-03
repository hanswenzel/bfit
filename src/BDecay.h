#ifndef BDECAY_H
#define BDECAY_H
#include <iostream>
#include "TObject.h"
#include "TNamed.h"
#include "Measurement.h"
#include <algorithm>
#include <vector>
#include <functional>
#include "bexclfit_fcts.h"
#include <TROOT.h>

 class BDecay : public TNamed 
{
 private:

  std::vector<Measurement>  sample;
  // cuts applied on sample
  double min_mass,max_mass,chiprob,chi2d,pt,ptk;
  int bitcode;
 
  //
  // The following class is used in the Prune method to apply more stringent mass requirements
  //
  class cut      // ! don't stream                         
    {
    public:
      cut(double min,double max);
      bool operator() (Measurement m)
	{
	  return !(m.Getmass()>min_mass&&m.Getmass()<max_mass);
	}
    private:
      double min_mass,max_mass;
    };
  //
  // The following class is used to select sideband
  //
  class sbcut    // ! don't stream
    {
    public:
      sbcut(double min,double max,double masswin);    
      bool operator() (Measurement m)
	{
	  return !((m.Getmass()<(min_mass+sbmasswin))||(m.Getmass()>(max_mass-sbmasswin)));
	}
    private:
      double min_mass,max_mass,sbmasswin;
    };

 public:

  BDecay():TNamed(){};
  BDecay(Text_t *name,const Text_t *title, const std::vector<Measurement>& vec):TNamed(name,title),sample(vec){};
  ~BDecay();
  void Prune(double mmin,double mmax);
  void Sort(){sort(sample.begin(),sample.end());};
  int Capacity(){return  sample.capacity();}
  int Size(){return  sample.size(); }
  std::vector<Measurement> *Sidebands(double min,double max,double masswin); // ! don't stream
  // iterators and iterator methods:
  std::vector<Measurement>::iterator enditer(){return sample.end();} 
  std::vector<Measurement>::iterator beginiter(){return sample.begin();}  
  // Accessor methods:
  const char * GetName() {return TNamed::fName;}
  // Accessor methods for cuts
  double Getmin_mass(){return min_mass;}
  double Getmax_mass(){return max_mass;}
  double Getchiprob(){return chiprob;}
  double Getchi2d(){return chi2d;}
  double Getpt_cut(){return pt;}
  double Getptk_cut(){return ptk;}
  int Getbitcode(){return bitcode;}
  void SetCuts(double, double, double, double, double, double); // required !
  
  ClassDef(BDecay,1) 
};
    
#endif   // BDECAY_H
    
