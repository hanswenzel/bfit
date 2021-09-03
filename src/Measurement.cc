#include "Measurement.h"
#include <iostream>

ClassImp(Measurement)
void   Measurement::Print() 
{ 
 std::cout << "ctau       :  " << ctau      << std::endl;
 std::cout << "sigma ctau :  " << sigmactau << std::endl; 
 std::cout << "mass       :  " << mass      << std::endl;
 std::cout << "sigmamass  :  " << sigmamass << std::endl;  
 std::cout << "iso        :  " << iso << std::endl;
 std::cout << "bpt        :  " << bpt << std::endl;
 std::cout << "kpt        :  " << kpt << std::endl;
}
 
Measurement::Measurement(const Measurement&  Measurement_in) // Copy Constructor 
{
  ctau      = Measurement_in.ctau;
  sigmactau = Measurement_in.sigmactau;
  mass      = Measurement_in.mass;
  sigmamass = Measurement_in.sigmamass;
  iso       = Measurement_in.iso;
  bpt       = Measurement_in.bpt;
  kpt       = Measurement_in.kpt;
}
bool Measurement::operator<(const Measurement&  Measurement_in)  const
{
  return  (mass<Measurement_in.mass);
}

