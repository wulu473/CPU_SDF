//Limits for small values
const static double areaLimit = 1.0E-8;
const static double zeroOffSet = 1.0E-6;

//PI
const static double PI = 3.14159265359;

//Replace unset inf values with the maxDistance of the apprpriate sign
void overWriteInfinityKernel(double* sdf, double maxValue, int length){

  for(int idx = 0; idx < length; idx++){
    float dist = sdf[idx].distance[1];
    if(isinf(dist)){
      if(dist > 0){
	sdf[idx].distance[1] = maxValue;
      }else{
	sdf[idx].distance[1] = -maxValue; 
      }
    }
  }
}

//Overwrite values that have amagnitude below some threshold
void overWriteSmallSDFKernel(double* sdf, double minValue, int length){

  for(int idx = 0; idx < length; idx++){
    if(fabs(sdf[idx]) < minValue){
      sdf[idx] = -minValue;
    }
  }
}

void fillUnsetValues(double* sdf, int width, int height){

//TODO
//Travel along rows and the columns from left -> right and then right -> left
//Once you see a transition from positive to negative between element n and n+1, set signs of following values to negative, until there is a transition to positive, in which case switch to positive sign

}



void getBoundingDimensions(double* rect_x, double* rect_y, double xMin, double yMin, double xMax, double yMax, double* boundMin, double* boundMax, unsigned int vertices){

  boundMin[0] = std::numeric_limits<double>::infinity();
  boundMin[1] = std::numeric_limits<double>::infinity();

  boundMax[0] = -std::numeric_limits<double>::infinity();
  boundMax[1] = -std::numeric_limits<double>::infinity();

  for(int i = 0; i < vertices; i++){
    boundMin[0] = min(boundMin[0], rect_x[i]);
    boundMin[1] = min(boundMin[1], rect_y[i]);

    boundMax[0] = max(boundMax[0], rect_x[i]);
    boundMax[1] = max(boundMax[1], rect_y[i]);
  }
}

void getRectangleCoords(double* edge_x, double* edge_y, double* rect_x, double* rect_y, double nx, double ny, double polygonExtent){

  //First two coordinates of the rectangle are just the edge co-ordinates
  rect_x[0] = edge_x[0]; 
  rect_y[0] = edge_y[0];
  rect_x[1] = edge_x[1]; 
  rect_y[1] = edge_y[1];

  //To get to 3rd and 4th co-ordinates of the rectangle, travel along the unit normal
  rect_x[2] = rect_x[1] + polygonExtent * nx; 
  rect_y[2] = rect_y[1] + polygonExtent * ny;
  rect_x[3] = rect_x[0] + polygonExtent * nx; 
  rect_y[3] = rect_y[0] + polygonExtent * ny;
}

void getTriangleCoords(double vertex_x, double vertex_y, double* trng_x, double* trng_y, double n0_x, double n0_y, double n1_x, double n1_y, double polygonExtent){

  trng_x[0] = vertex_x;
  trng_y[0] = vertex_y;

  trng_x[1] = trng_x[0] + polygonExtent * n0_x; 
  trng_y[1] = trng_y[0] + polygonExtent * n0_y;

  trng_x[2] = trng_x[0] + polygonExtent * n1_x; 
  trng_y[2] = trng_y[0] + polygonExtent * n1_y;
}

/**
 *	Populate local cells with signed distance to geometry using the CSC algorithm
 *	@param sdf pre-allcoated row-major pointer to SDF data at cell centres
 *	@param xMin value at the leftmost edge of the whole bounding volume
 *	@param yMin value at the bottommost edge of the whole bounding volume
 *	@param dx cell length in x-dimension
 *	@param dy cell length in y-dimension
 *	@param left offset of leftmost cell wrt bounding volume
 *	@param right offset of rightmost+1 cell wrt bounding volume
 *	@param bottom offset of bottommost cell wrt bounding volume
 *	@param top offset of topmost+1 cell wrt bounding volume
 *	@param vertices list of vertices in geometry [x0, y0, x1, y1, ... x(n-1), y(n-1)]
 *	@param numVertices number of vertices in geometry
 *	@param maxDist maximum distance of SDF
 */
void getSDF(double* sdf, double xMin, double yMin, double dx, double dy, int left, int right, int bottom, int top, const double* vertices, int numVertices, double maxDist){

  int width = right - left;
  int height = top - bottom;

  //Number of cells in local domain
  int length = width * height;

  //Maximum limits of the local domain
  double xMax = xMin + width * dx;
  double yMax = yMin + height * dy;

  //Initialise all values in sdf to inf
  for(int i = 0; i < length; i++){
    sdf[i] = std::numeric_limits<double>::infinity();
  }

  //Coordinates of a vertex
  double vertex_x;
  double vertex_y;

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
  double* outwardNormals = malloc((2 * numVertices) * sizeof(double));

  //Extrusion bounding volume
  double boundMin[2];
  double boundMax[2];

  //Dimensions of extrusion bounding volume
  int width;
  int height;

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

    //Calculate inward unit normal for edge
    double in_norm_x = -out_norm_x; 
    double in_norm_y = -out_norm_y;

    //Get positive edge extrusion coordinates
    getRectangleCoords(edge_x, edge_y, rect_x, rect_y, out_norm_sin, out_norm_cos, polygonExtent);

    //Get the interection of the extrusion boiunding volume and the domain
    getBoundingDimensions(rect_x, rect_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 4);

    boundMin[0] = floor((max(xMin, boundMin[0])-xMin)/dx); 
    boundMin[1] = floor((max(yMin, boundMin[1])-yMin)/dy); 

    boundMax[0] = ceil((min(xMax, boundMax[0])-xMin)/dx); 
    boundMax[1] = ceil((min(yMax, boundMax[1])-yMin)/dy);  

    width = boundMax[0] - boundMin[0];
    height = boundMax[1] - boundMin[1];

    //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
    if((width > 0) && (height > 0)){
      //Update positive distances for grid points in edge rectangle
      updateRectangle(true, sdf, rect_x, rect_y);
    }

    //Get negative edge extrusion coordinates
    getRectangleCoords(edge_x, edge_y, rect_x, rect_y, in_norm_sin, in_norm_cos, polygonExtent);

    //Get the interection of the extrusion boiunding volume and the domain
    getBoundingDimensions(rect_x, rect_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 4);

    boundMin[0] = floor((max(xMin, boundMin[0])-xMin)/dx); 
    boundMin[1] = floor((max(yMin, boundMin[1])-yMin)/dy); 

    boundMax[0] = ceil((min(xMax, boundMax[0])-xMin)/dx); 
    boundMax[1] = ceil((min(yMax, boundMax[1])-yMin)/dy);  

    width = boundMax[0] - boundMin[0];
    height = boundMax[1] - boundMin[1];

    //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
    if((width > 0) && (height > 0)){
      //Update negative distances for grid points in edge rectangle
      updateRectangle(false, sdf, rect_x, rect_y);
    }
  }

  //sdf has distance information from all edges

  //-------------------- VERTEX --------------------

  //Vector to previous vertex
  double toPrev[2];
  //Vector to next vertex
  double toNext[2];

  //Loop over vertices, constructing +ve and -ve trangles and scan convert the grid points
  for(int i = 0; i < (numVertices); i++){

    //Indices of adjacent vertices (vertices are globally ordered clockwise)
    int previousVertex = i - 1;
    if(i == 0){
      previousVertex = numVertices - 1;
    }

    int nextVertex = i + 1;
    if(i == (numVertices - 1)){
      previousVertex = 0;
    }

    //Construct adjacent edge vectors to calculate the vertex angle 
    toPrev[0] = vertices[2*previous_vertex_index] - vertices[2*i];
    toPrev[1] = vertices[2*previous_vertex_index+1] - vertices[2*i+1];

    toNext[0] = vertices[2*next_vertex_index] - vertices[2*i];
    toNext[1] = vertices[2*next_vertex_index+1] - vertices_y[2*i+1];

    //Calculate vertex angle. atan2, cross and dot products are combined to ensure we get a sense of the sign of the vertex angle
    double vertexAngleDeg = atan2(((toNext[0] * toPrev[1]) - (toPrev[0] * toNext[1])), (toPrev[0]*toNext[0] + toPrev[1]*toNext[1])) * 180.0 / PI;

    //atan2 returns values in range [-pi, pi] with a discontinuity at the transition. Adding 360 degrees to negative angles gives values of the whole  unit circle. 
    if(vertexAngleDeg < 0){
      vertexAngleDeg = vertexAngleDeg + 360.0;
    }

    //Vertices with angles > 180 degrees are convex and require positive extrusions. Angles < 180 degrees denote concave vertices requiring negative extrusions. Angles of 180 degrees are flat and need no extrusions
    if(vertexAngleDeg > 180.0){
      //Get positive distance triangle coordinates
      getTriangleCoords(vertices[2*i], vertices[2*i+1], trng_x, trng_y, outwardNormals[2*previousVertex], outwardNormals[2*previousVertex+1], outwardNormals[2*i], outwardNormals[2*i+1], maxDist);

       //Get the interection of the extrusion boiunding volume and the domain
      getBoundingDimensions(trng_x, trng_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 3);

      boundMin[0] = floor((max(xMin, boundMin[0])-xMin)/dx); 
      boundMin[1] = floor((max(yMin, boundMin[1])-yMin)/dy); 

      boundMax[0] = ceil((min(xMax, boundMax[0])-xMin)/dx); 
      boundMax[1] = ceil((min(yMax, boundMax[1])-yMin)/dy);  

      width = boundMax[0] - boundMin[0];
      height = boundMax[1] - boundMin[1];

      //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
      if((width > 0) && (height > 0)){
        //Update positive distances for grid points in vertex rectangle
        updateTriangle(true, sdf, trng_x, trng_y);
      }

    }else if(vertexAngleDeg < 180.0){
      //Get negative distance triangle coordinates
      getTriangleCoords(vertices[2*i], vertices[2*i+1], trng_x, trng_y, -outwardNormals[2*previousVertex], -outwardNormals[2*previousVertex+1], -outwardNormals[2*i], -outwardNormals[2*i+1], maxDist);
     
      //Get the interection of the extrusion boiunding volume and the domain
      getBoundingDimensions(trng_x, trng_y, xMin, yMin, xMax, yMax, boundMin, boundMax, 3);

      boundMin[0] = floor((max(xMin, boundMin[0])-xMin)/dx); 
      boundMin[1] = floor((max(yMin, boundMin[1])-yMin)/dy); 

      boundMax[0] = ceil((min(xMax, boundMax[0])-xMin)/dx); 
      boundMax[1] = ceil((min(yMax, boundMax[1])-yMin)/dy);  

      width = boundMax[0] - boundMin[0];
      height = boundMax[1] - boundMin[1];

      //If there is overlap between local and extrusion bounding volumes, do distance calculation for cells in that intersection
      if((width > 0) && (height > 0)){
        //Update positive distances for grid points in vertex rectangle
        updateTriangle(false, sdf, trng_x, trng_y);
      }
    }
  }

  //sdf has distance values from all edge and vertex extrusions

  //Sweep along rows and columns to fill in unset values
  fillUnsetValues(sdf, width, height);
  

  //Overwrite infinity with maxDist of correct signe
  overWriteInfinityKernel(sdf, maxValue, length);
    
  //Optional
  overWriteSmallSDFKernel(sdf, minValue, length);

  free(outwardNormals);
}
