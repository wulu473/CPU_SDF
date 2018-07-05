#include <stdlib.h>
#include <limits>
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include "2DSDF.hpp"
#include "SignedDistance.hpp"

int main(int argc, char** argv){

  int width = 20;
  int height = 10;

  int length = width * height;

  double xMin = -10.0;
  double yMin =  -7.0;

  double xMax =  10.0;
  double yMax =   3.0;

  int left = 0;
  int right = 0;
  int bottom = 0;
  int top = 0;

  double dx = (xMax - xMin) / (double)width;
  double dy = (yMax - yMin) / (double)height;

  double maxDistance = 5.0;

  double vertices[6] = {-1,0,
                         0,1,
                         1,0};

  int numVertices = 3;

  double* sdf = (double*)malloc(length * sizeof(double));
  getSDF(sdf, xMin, yMin, dx, dy, left, right, bottom, top, vertices, numVertices, maxDistance);

  for(int y = 0; y < width; y++){
    for(int x = 0; x < width; x++){
    
      printf("%f\t", sdf[y * width + x]);

    }
    printf("\n");
  }

  return 0;
}

