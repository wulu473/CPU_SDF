#include <stdio.h>
#include <iostream>
#include <fstream>
#include "SignedDistance.hpp"

int main(int argc, char** argv){

  int width = 200;
  int height = 100;

  int length = width * height;

  double xMin = -10.0;
  double yMin =  -5.0;

  double xMax =  10.0;
  double yMax =   5.0;

  double dx = (xMax - xMin) / (double)width;
  double dy = (yMax - yMin) / (double)height;

  double maxDistance = 0.5;

  double vertices[8] = {-4,-2,
                         0,3,
                         4,-2,
                         0, -1};

  int numVertices = 4;

  double* sdf = (double*)malloc(length * sizeof(double));
  getSDF(sdf, xMin, yMin, dx, dy, width, height, vertices, numVertices, maxDistance);

  //Output data file
  std::ofstream myfile;
  myfile.open ("sdf.dat");

  for(int y = 0; y < height; y++){
    for(int x = 0; x < width; x++){

      myfile << sdf[y * width + x] << " ";

    }
    myfile << "\n";
  }

  myfile.close();

  return 0;
}

