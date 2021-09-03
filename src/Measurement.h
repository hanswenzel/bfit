#ifndef MEASUREMENT_H
#define MEASUREMENT_H
#include "TObject.h"
class Measurement 
{ 
 private:
      double ctau;
      double sigmactau;
      double mass;
      double sigmamass;  
      double iso;
      double bpt;
      double kpt;
 public:
      // Default Constructor set everything to 0 
      Measurement(){ctau=0.0;sigmactau=0.0;mass=0.0;sigmamass=0.0;iso=0.0;bpt=0.0;kpt=0.0;};
      // Constructor initializing everything  
      Measurement(double l, double s, double m, double sm, double i, double b, double k): ctau(l),sigmactau(s),mass(m),sigmamass(sm),iso(i),bpt(b),kpt(k) {}; 

      // Constructor initializing all 
      void   Print(); 
      double Getctau(){return ctau;}
      double Getsigmactau(){return sigmactau;}
      double Getmass(){return mass;}
      double Getsigmamass(){return sigmamass;} 
      double Getiso(){return iso;} 
      double Getpt(){return bpt;}
      double Getkpt(){return kpt;}
      bool operator<(const Measurement&) const; // test m1<m2
      Measurement(const Measurement&); // Copy Constructor
      ClassDef(Measurement,1)
};

#endif   // MEASUREMENT_H
