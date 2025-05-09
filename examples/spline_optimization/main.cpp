#include <array>
#include <iostream>
#include "irl/splines/Spline.h"
#include <math.h>
#include <nlopt.hpp>



int main() {
  int Nd = 12;
  std::cout << "Start\n";
  nlopt::opt opt(nlopt::LN_COBYLA,Nd);
  
  // Define Rectangle for Clipping
  std::vector<std::vector<double>> square = {
    {0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}, {0.0, 0.0}};

  // Pick Interpolating Points and Tangents (Vibes)
  std::vector<std::vector<double>> InterpolatingPoints = {
    {-1.0, 0.0}, {0.5, -1.0}, {1.5, 0.5},
    {0.75, 0.6}, {0.25, 0.5}, {-1.0, 0.0}};
  std::vector<std::vector<double>> Tangents = {
    {-0.555, -0.832}, {1.0, 0},        {-0.196, 0.981},
    {-0.707, -0.707}, {-0.707, 0.707}, {-0.555, -0.832}};
    
  // Circle
  double R = 0.999999;
  std::vector<std::vector<double>> Circle =
  {{R,0.0},{0.0,R},{-R,0.0},{0.0,-R},{R,0.0}};
  std::vector<std::vector<double>> CircleT =
  {{0.0,1.0},{-1.0,0.0},{0.0,-1.0},{1.0,0.0},{0.0,1.0}};


  // InterpolatingPoints = Circle;
  // Tangents = CircleT;
  // Interpolate
  Spline s = Spline<double>::LocalRQuadInterp(InterpolatingPoints, Tangents);

  return 0;
}