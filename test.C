#include <cmath>
#include <cassert>
#include <limits>

void getSDF(double* sdf, double xMin, double yMin,
    double dx, double dy, int left, int right, int bottom, int top,
    const double* vertices, int numVertices, double maxDist);

int main(int, char **)
{
  double triangle[6] = {-1,0,
                         0,1,
                         1,0};

  // Test 1
  // Check if a point on an edge of a triangle has distanc e0
  double sdf = std::numeric_limits<double>::signaling_NaN();
  getSDF(&sdf,-0.5,-0.5,1.0,1.0,0,0,1,1,triangle,3,100);
  assert(fabs(sdf) < 1e-10);

  return 0;
}

