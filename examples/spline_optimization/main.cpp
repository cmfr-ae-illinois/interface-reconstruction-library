#include <array>
#include <iostream>
#include "irl/splines/Spline.h"
#include <math.h>
#include <nlopt.hpp>



int main() {
  int Nd = 12;
  std::cout << "Start\n";
  nlopt::opt opt(nlopt::LN_COBYLA,Nd);
  
  

  return 0;
}