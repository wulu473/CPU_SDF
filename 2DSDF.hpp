//Half plane test. Point is within rectangular extrusion if it is on the same side of all its edges
bool pointInTriangle(double x, double y, double* trng_x, double* trng_y, double* normalA, double* normalB, double* normalC, double limit){

  if((normalA[0] * (trng_x[0]-x) + (normalA[1]) * (trng_y[0]-y))>(limit)){
    return false;
  }

  if((normalB[0] * (trng_x[1]-x) + (normalB[1]) * (trng_y[1]-y))>(limit)){
    return false;
  }

  if((normalC[0] * (trng_x[0]-x) + (normalC[1]) * (trng_y[0]-y))>(limit)){
    return false;
  }

  return true;
}

//Half plane test. Point is within rectangular extrusion if it is on the same side of all its edges
bool pointInRectangle(double x, double y, double* rect_x, double* rect_y, double* normalA, double* normalB, double* normalC, double* normalD, double limit){

  if((normalA[0] * (rect_x[0]-x) + (normalA[1]) * (rect_y[0]-y))>(limit)){
    return false;
  }

  if((normalB[0] * (rect_x[1]-x) + (normalB[1]) * (rect_y[1]-y))>(limit)){
    return false;
  }

  if((normalC[0] * (rect_x[2]-x) + (normalC[1]) * (rect_y[2]-y))>(limit)){
    return false;
  }

  if((normalD[0] * (rect_x[0]-x) + (normalD[1]) * (rect_y[0]-y))>(limit)){
    return false;
  }

  return true;
}

void getLineNormal(double* norm, double aX, double aY, double bX, double bY){
  norm[0] = -(bY - aY);
  norm[1] = bX - aX;

  double mag = std::sqrt((norm[0]*norm[0])+(norm[1]*norm[1]));
  norm[0] /= mag;
  norm[1] /= mag;
}

//Loop trhough cells in rectangular region of the mesh and determine if cell centres are within an edge extrusion. If so, find the distance to the edge and record the minimum magnitude distance
void setEdgeSDF(bool positive, double* sdf, double* rect_x, double* rect_y, double xMin, double yMin, int startX, int startY, int endX, int endY, int width, double dx, double dy, double maxDistance, double limit){

  //Inward pointing normals of the extrusion
  double normA[2];
  double normB[2];
  double normC[2];
  double normD[2];

  getLineNormal(normA, rect_x[0], rect_y[0], rect_x[1], rect_y[1]);
  getLineNormal(normB, rect_x[1], rect_y[1], rect_x[2], rect_y[2]);
  getLineNormal(normC, rect_x[2], rect_y[2], rect_x[3], rect_y[3]);
  getLineNormal(normD, rect_x[3], rect_y[3], rect_x[0], rect_y[0]);

  //Coordinates of current cell
  double x_;
  double y_;

  //For all cells in region
  for(int y = startY; y <= endY; y++){
    for(int x = startX; x <= endX; x++){

      x_ = xMin + x * dx + 0.5 * dx;
      y_ = yMin + y * dy + 0.5 * dy;

      //If cell is in extrusion
      if(pointInRectangle(x_, y_, rect_x, rect_y, normA, normB, normC, normD, limit)){
	//Find distance to edge
	double dist = (normA[0] * (rect_x[0]-x_) + (normA[1]) * (rect_y[0]-y_));
	if(positive){
	  dist = -dist;
	}
	//If distance has small enough magnitude
	if(std::abs(dist) < maxDistance){
	  //If distance magnitude is smaller than previous value
	  if(std::abs(sdf[y * width + x]) > std::abs(dist)){
	    //Write distance to memory with appropriate sign
	    sdf[y * width + x] = dist;
	  }
	}
      }
    }
  } 
}

//Loop trhough cells in rectangular region of the mesh and determine if cell centres are within a vertex extrusion. If so, find the distance to the vertex and record the minimum magnitude distance
void setVertexSDF(bool positive, double* sdf, double* trng_x, double* trng_y, double xMin, double yMin, int startX, int startY, int endX, int endY, int width, double dx, double dy, double maxDistance, double limit){

  //Inward pointing normals of the extrusion
  double normA[2];
  double normB[2];
  double normC[2];

  //TODO 
  //Check edge order
  getLineNormal(normA, trng_x[1], trng_y[1],trng_x[0], trng_y[0]);
  getLineNormal(normB, trng_x[2], trng_y[2],trng_x[1], trng_y[1]);
  getLineNormal(normC, trng_x[0], trng_y[0],trng_x[2], trng_y[2]);

  //Coordinates of current cell
  double x_;
  double y_;
 
  //Sign of distance 
  int sign;

  if(positive){
    sign = 1;
  }else{
    sign = -1;
  }

  //For all cells in region
  for(int y = startY; y <= endY; y++){
    for(int x = startX; x <= endX; x++){

      x_ = xMin + x * dx + 0.5 * dx;
      y_ = yMin + y * dy + 0.5 * dy;

      //If cell is in extrusion
      if(pointInTriangle(x_, y_, trng_x, trng_y, normA, normB, normC, limit)){
	//Find distance to vertex
	double dist = std::sqrt((x_ - trng_x[0]) * (x_ - trng_x[0]) + (y_ - trng_y[0]) * (y_ - trng_y[0]));
	//If distance has small enough magnitude
	if(dist < maxDistance){
	  //If distance magnitude is smaller than previous value
	  if(std::abs(sdf[y * width + x]) > dist){
	    //Write distance to memory with appropriate sign
	    sdf[y * width + x] = sign * dist;
	  }
	}
      } 
    }
  }
}
