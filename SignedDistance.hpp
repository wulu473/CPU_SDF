#include <stdlib.h>
#include <limits>
#include <cmath>
#include <algorithm>
#include "2DSDF.hpp"

//PI
const static double PI = 3.14159265359;

//Replace unset inf values with the maxDistance of the apprpriate sign
void overWriteInfinityKernel(double* sdf, double maxValue, int length){

  for(int idx = 0; idx < length; idx++){
    float dist = sdf[idx];
    if(std::isinf(dist)){
      if(dist > 0){
	sdf[idx] = maxValue;
      }else{
	sdf[idx] = -maxValue; 
      }
    }
  }
}

//Overwrite values that have amagnitude below some threshold
void overWriteSmallSDFKernel(double* sdf, double minValue, int length){

  for(int idx = 0; idx < length; idx++){
    if(std::abs(sdf[idx]) < minValue){
      sdf[idx] = -minValue;
    }
  }
}

//Sweep along rows and columns to fill in unset values in the interior of the geometry
void fillUnsetValues(double* sdf, int width, int height){

  double previous = 1.0/0.0;
  double value;
  bool flip;

  for(int y = 0; y < height; y++){
    //Left to right
    previous = 1.0/0.0;
    flip = false;
    for(int x = 0; x < width; x++){
      value = sdf[y * width + x];

      if((previous < 0.0) && (std::isinf(value)) && !(std::isinf(previous))){
	flip = true;
      }

      if((std::isinf(previous) && !(std::isinf(value)) && (value < 0.0))){
	flip = false;
      }

      previous = value;

      if((flip)&&(value > 0.0)){
	//Fill with -inf
	sdf[y * width + x] = -1.0/0.0;
      }
    }

    //Right to left
    previous = 1.0/0.0;
    flip = false;
    for(int x = (width-1); x >= 0; x--){
      value = sdf[y * width + x];

      if((previous < 0.0) && (std::isinf(value)) && !(std::isinf(previous))){
	flip = true;
      }

      if((std::isinf(previous) && !(std::isinf(value)) && (value < 0.0))){
	flip = false;
      }

      previous = value;

      if((flip)&&(value > 0.0)){
	//Fill with -inf
	sdf[y * width + x] = -1.0/0.0;
      }
    }
  }

  for(int x = 0; x < width; x++){
    //Bottom to top
    previous = 1.0/0.0;
    flip = false;
    for(int y = 0; y < height; y++){
      value = sdf[y * width + x];

      if((previous < 0.0) && (std::isinf(value)) && !(std::isinf(previous))){
	flip = true;
      }

      if((std::isinf(previous) && !(std::isinf(value)) && (value < 0.0))){
	flip = false;
      }

      previous = value;

      if((flip)&&(value > 0.0)){
	//Fill with -inf
	sdf[y * width + x] = -1.0/0.0;
      }
    }

    //Top to bottom
    previous = 1.0/0.0;
    flip = false;
    for(int y = (height-1); y >= 0; y--){
      value = sdf[y * width + x];

      if((previous < 0.0) && (std::isinf(value)) && !(std::isinf(previous))){
	flip = true;
      }

      if((std::isinf(previous) && !(std::isinf(value)) && (value < 0.0))){
	flip = false;
      }

      previous = value;

      if((flip)&&(value > 0.0)){
	//Fill with -inf
	sdf[y * width + x] = -1.0/0.0;
      }
    }
  }
}

//Find the axis aligned bounding box around an extrusion and its intersection with the domain
void getBoundingDimensions(double* poly_x, double* poly_y, double xMin, double yMin, double xMax, double yMax, double* boundMin, double* boundMax, unsigned int vertices, double maxDistance){

  boundMin[0] = std::numeric_limits<double>::infinity();
  boundMin[1] = std::numeric_limits<double>::infinity();

  boundMax[0] = -std::numeric_limits<double>::infinity();
  boundMax[1] = -std::numeric_limits<double>::infinity();

  for(unsigned int i = 0; i < vertices; i++){
    boundMin[0] = std::min(boundMin[0], poly_x[i]);
    boundMin[1] = std::min(boundMin[1], poly_y[i]);

    boundMax[0] = std::max(boundMax[0], poly_x[i]);
    boundMax[1] = std::max(boundMax[1], poly_y[i]);
  }

  //Limit a vertex bouding box by finding its intersection with the bounding box of a disc centred at the vertex
  if(vertices == 3){
    boundMin[0] = std::max(boundMin[0], poly_x[0] - maxDistance);
    boundMin[1] = std::max(boundMin[1], poly_y[0] - maxDistance);

    boundMax[0] = std::min(boundMax[0], poly_x[0] + maxDistance);
    boundMax[1] = std::min(boundMax[1], poly_y[0] + maxDistance);
  }

  boundMin[0] = std::max(boundMin[0], xMin);
  boundMin[1] = std::max(boundMin[1], yMin);

  boundMax[0] = std::min(boundMax[0], xMax);
  boundMax[1] = std::min(boundMax[1], yMax);

}

//Find the extrusion of an edge. This is constructed by projecting the edge vertices in the normal direction. The vertices of the rectangle must be in a CCW order to produce inward pointing normals in the half-plane test later.
void getRectangleCoords(bool positive, double* edge_x, double* edge_y, double* rect_x, double* rect_y, double nx, double ny, double polygonExtent){

  //First two coordinates of the rectangle are just the edge co-ordinates
  if(positive){
    rect_x[0] = edge_x[0]; 
    rect_y[0] = edge_y[0];
    rect_x[1] = edge_x[1]; 
    rect_y[1] = edge_y[1];
  }else{
    rect_x[0] = edge_x[1]; 
    rect_y[0] = edge_y[1];
    rect_x[1] = edge_x[0]; 
    rect_y[1] = edge_y[0];
  }

  //To get to 3rd and 4th co-ordinates of the rectangle, travel along the unit normal
  rect_x[2] = rect_x[1] + polygonExtent * nx; 
  rect_y[2] = rect_y[1] + polygonExtent * ny;
  rect_x[3] = rect_x[0] + polygonExtent * nx; 
  rect_y[3] = rect_y[0] + polygonExtent * ny;
}

//Find the extrusion of a vertex. This is constructed by projecting the vertex in the normal directions of the edges that meet at the vertex. The vertices of the triangle must be in a CCW order to produce inward pointing normals in the half-plane test later
void getTriangleCoords(double vertex_x, double vertex_y, double* trng_x, double* trng_y, double n0_x, double n0_y, double n1_x, double n1_y, double maxDistance){

  //Find the angle between the normals
  double angle = atan2(((n1_x * n0_y) - (n0_x * n1_y)), (n0_x*n1_x + n0_y*n1_y));

  //Side of triangle with height maxDistance
  double polygonExtent = maxDistance / std::cos(angle/2.0); 

  //The first element of the triangle is the vertex itself
  trng_x[0] = vertex_x;
  trng_y[0] = vertex_y;

  trng_x[1] = trng_x[0] + polygonExtent * n1_x; 
  trng_y[1] = trng_y[0] + polygonExtent * n1_y;

  trng_x[2] = trng_x[0] + polygonExtent * n0_x; 
  trng_y[2] = trng_y[0] + polygonExtent * n0_y;
}

/**
 *	Populate local cells with signed distance to geometry using the CSC algorithm
 *	@param sdf pre-allcoated row-major pointer to SDF data at cell centres
 *	@param xMin value at the leftmost edge of the local region
 *	@param yMin value at the bottommost edge of the local region
 *	@param dx cell length in x-dimension
 *	@param dy cell length in y-dimension
 *	@param width number of cells in x-dimension
 *	@param height number of cells in y-dimension
 *	@param vertices list of vertices in geometry [x0, y0, x1, y1, ... x(n-1), y(n-1)]
 *	@param numVertices number of vertices in geometry
 *	@param maxDist maximum distance of SDF
 */
void getSDF(double* sdf, double xMin, double yMin, double dx, double dy, int width, int height, const double* vertices, int numVertices, double maxDistance){

  //Number of cells in local domain
  int length = width * height;

  //Local incideces for the limits of extrusion bounding volumes
  int startX;
  int startY;

  int endX;
  int endY;

  //Maximum limits of the local domain
  double xMax = xMin + width * dx;
  double yMax = yMin + height * dy;

  //Initialise all values in sdf to inf
  for(int i = 0; i < length; i++){
    sdf[i] = std::numeric_limits<double>::infinity();
  }

  //Epsilon for point-in-polygon test
  double limit = 0.01 * std::min(dx, dy);

  //Coordinates of vertices forming an edge
  double edge_x[2];
  double edge_y[2];

  //Coordinates of edge extrusion corners
  double rect_x[4];
  double rect_y[4];

  //Coordinates of vertex extrusion corners
  double trng_x[3];
  double trng_y[3];

  //Outward pointing edge normals
  double* outwardNormals = (double*)malloc((2 * numVertices) * sizeof(double));

  //Extrusion bounding volume
  double boundMin[2];
  double boundMax[2];

  //Dimensions of extrusion bounding volume
  int b_width;
  int b_height;

  //sdf is infinity everywhere

  //-------------------- EDGE --------------------

  //Loop over edges, constructing +ve and -ve rectangles and scan convert the grid points
  for(int i = 0; i < (numVertices * 2); i+=2){

    edge_x[0] = vertices[i]; 
    edge_y[0] = vertices[i+1]; 

    //Final edge is last vertex connected to first 
    if(i < ((numVertices * 2) - 2)){
      edge_x[1] = vertices[i+2];
      edge_y[1] = vertices[i+3];
    }else{
      edge_x[1] = vertices[0];
      edge_y[1] = vertices[1];
    }

    //Calculate outward unit normal for edge (edge points assumed to be in clockwise order)
    double edge_dx = edge_x[1] - edge_x[0]; 
    double edge_dy = edge_y[1] - edge_y[0];

    double out_norm_mag = sqrt(edge_dx*edge_dx + edge_dy*edge_dy);
    double out_norm_x = (-edge_dy)/out_norm_mag; 
    double out_norm_y = (edge_dx)/out_norm_mag;

    //Store edge normals for vertex extrusions
    outwardNormals[i]   = out_norm_x;
    outwardNormals[i+1] = out_norm_y;

    //Get positive edge extrusion coordinates
    getRectangleCoords(true, edge_x, edge_y, rect_x, rect_y, out_norm_x, out_norm_y, maxDistance);

    //Get the interection of the extrusion bounding volume and the domain
    getBoundingDimensions(rect_x, rect_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 4, maxDistance);

    //Find the indices of cells at the edges of the bounding volume intersection
    startX = std::floor((boundMin[0]-xMin)/dx); 
    startY = std::floor((boundMin[1]-yMin)/dy); 

    endX = std::ceil((boundMax[0]-xMin)/dx); 
    endY = std::ceil((boundMax[1]-yMin)/dy);  

    b_width = endX - startX;
    b_height = endY - startY;

    //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
    if((b_width > 0) && (b_height > 0)){
      //Update positive distances for grid points in edge rectangle
      setEdgeSDF(true, sdf, rect_x, rect_y, xMin, yMin, startX, startY, endX, endY, width, dx, dy, maxDistance, limit);
    }

    //Get negative edge extrusion coordinates
    getRectangleCoords(false, edge_x, edge_y, rect_x, rect_y, -out_norm_x, -out_norm_y, maxDistance);

    //Get the interection of the extrusion bounding volume and the domain
    getBoundingDimensions(rect_x, rect_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 4, maxDistance);

    //Find the indices of cells at the edges of the bounding volume intersection
    startX = std::floor((boundMin[0]-xMin)/dx); 
    startY = std::floor((boundMin[1]-yMin)/dy); 

    endX = std::ceil((boundMax[0]-xMin)/dx); 
    endY = std::ceil((boundMax[1]-yMin)/dy);  

    b_width = endX - startX;
    b_height = endY - startY;

    //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
    if((b_width > 0) && (b_height > 0)){
      //Update negative distances for grid points in edge rectangle
      setEdgeSDF(false, sdf, rect_x, rect_y, xMin, yMin, startX, startY, endX, endY, width, dx, dy, maxDistance, limit);
    }
  }

  //sdf has distance information from all edges

  //-------------------- VERTEX --------------------

  //Vector to previous vertex
  double toPrev[2];
  //Vector to next vertex
  double toNext[2];

  //Loop over vertices, constructing +ve and -ve trangles and scan convert the grid points
  for(int i = 0; i < numVertices; i++){

    //Indices of adjacent vertices (vertices are globally ordered clockwise)
    int previousVertex = i - 1;
    if(i == 0){
      previousVertex = numVertices - 1;
    }

    int nextVertex = i + 1;
    if(i == (numVertices - 1)){
      nextVertex = 0;
    }

    //Construct adjacent edge vectors to calculate the vertex angle 
    toPrev[0] = vertices[2*previousVertex] - vertices[2*i];
    toPrev[1] = vertices[2*previousVertex+1] - vertices[2*i+1];

    toNext[0] = vertices[2*nextVertex] - vertices[2*i];
    toNext[1] = vertices[2*nextVertex+1] - vertices[2*i+1];

    //Calculate vertex angle. atan2, cross and dot products are combined to ensure we get a sense of the sign of the vertex angle
    double vertexAngleDeg = atan2(((toNext[0] * toPrev[1]) - (toPrev[0] * toNext[1])), (toPrev[0]*toNext[0] + toPrev[1]*toNext[1])) * 180.0 / PI;

    //atan2 returns values in range [-pi, pi] with a discontinuity at the transition. Adding 360 degrees to negative angles gives values of the whole  unit circle. 
    if(vertexAngleDeg < 0){
      vertexAngleDeg = vertexAngleDeg + 360.0;
    }

    //Vertices with angles > 180 degrees are convex and require positive extrusions. Angles < 180 degrees denote concave vertices requiring negative extrusions. Angles of 180 degrees are flat and need no extrusions
    if(vertexAngleDeg > 180.0){
      //Get positive distance triangle coordinates
      getTriangleCoords(vertices[2*i], vertices[2*i+1], trng_x, trng_y, outwardNormals[2*previousVertex], outwardNormals[2*previousVertex+1], outwardNormals[2*i], outwardNormals[2*i+1], maxDistance);

      //Get the interection of the extrusion bounding volume and the domain
      getBoundingDimensions(trng_x, trng_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 3, maxDistance);

      //Find the indices of cells at the edges of the bounding volume intersection
      startX = std::floor((boundMin[0]-xMin)/dx); 
      startY = std::floor((boundMin[1]-yMin)/dy); 

      endX = std::ceil((boundMax[0]-xMin)/dx); 
      endY = std::ceil((boundMax[1]-yMin)/dy);  

      b_width = endX - startX;
      b_height = endY - startY;

      //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
      if((b_width > 0) && (b_height > 0)){
	//Update positive distances for grid points in vertex rectangle
	setVertexSDF(true, sdf, trng_x, trng_y, xMin, yMin, startX, startY, endX, endY, width, dx, dy, maxDistance, limit);
      }

    }else if(vertexAngleDeg < 180.0){
      //Get negative distance triangle coordinates
      //previous->i->next is switched wrt positive normals
      getTriangleCoords(vertices[2*i], vertices[2*i+1], trng_x, trng_y, -outwardNormals[2*i], -outwardNormals[2*i+1], -outwardNormals[2*previousVertex], -outwardNormals[2*previousVertex+1], maxDistance);

      //Get the interection of the extrusion bounding volume and the domain
      getBoundingDimensions(trng_x, trng_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 3, maxDistance);

      //Find the indices of cells at the edges of the bounding volume intersection
      startX = std::floor((boundMin[0]-xMin)/dx); 
      startY = std::floor((boundMin[1]-yMin)/dy); 

      endX = std::ceil((boundMax[0]-xMin)/dx); 
      endY = std::ceil((boundMax[1]-yMin)/dy);  

      b_width = endX - startX;
      b_height = endY - startY;

      //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
      if((b_width > 0) && (b_height > 0)){
	//Update negative distances for grid points in vertex rectangle
	setVertexSDF(false, sdf, trng_x, trng_y, xMin, yMin, startX, startY, endX, endY, width, dx, dy, maxDistance, limit);
      }
    }
  }

  //sdf has distance values from all edge and vertex extrusions

  //Sweep along rows and columns to fill in unset values
  fillUnsetValues(sdf, width, height);  

  //Overwrite infinity with maxDist of correct sign
  overWriteInfinityKernel(sdf, maxDistance, length);

  //Optional
  //overWriteSmallSDFKernel(sdf, minValue, length);

  free(outwardNormals);
}
