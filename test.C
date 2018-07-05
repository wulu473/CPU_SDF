#define BOOST_TEST_MODULE Main
#define BOOST_TEST_DYN_LINK

#include <cmath>
#include <cassert>
#include <limits>
#include <boost/test/unit_test.hpp>

void getSDF(double* sdf, double xMin, double yMin,
    double dx, double dy, int left, int right, int bottom, int top,
    const double* vertices, int numVertices, double maxDist);

BOOST_AUTO_TEST_SUITE(SignedDistanceTest)

BOOST_AUTO_TEST_CASE(PolygonTest)
{
  double poly[10] = {0,0,
                     0,2,
                     1,2,
                     2,1,
                     1,0};

  auto distance = [poly](const double x, const double y)
  {
    double sdf = std::numeric_limits<double>::signaling_NaN();
    getSDF(&sdf,x-0.5,y-0.5,1.0,1.0,0,0,1,1,poly,5,100);
    return sdf;
  };

  // Point inside polygon OK
  BOOST_CHECK_GT(distance(0.5,1), 0);

  BOOST_CHECK_CLOSE(distance(0.5,1), 0.5, 1e-5);

  // Point inside polygon OK
  BOOST_CHECK_GT(distance(1.5,1), 0);

  // Shortest distance goes along line perp. to sloping side
  BOOST_CHECK_CLOSE(distance(1.5,1), 0.5 * sqrt(2*0.5*0.5), 1e-5);
  
  // Point outside polygon OK
  BOOST_CHECK_LT(distance(2,0), 0 );

  // Shortest distance goes along normal to sloping side
  BOOST_CHECK_CLOSE(distance(2,0), - 0.5 * sqrt(2*1*1), 1e-5 );
  
  // Point outside polygon OK
  BOOST_CHECK_LT(distance(0.5, 2.1), 0 );

  BOOST_CHECK_CLOSE(distance(0.5, 2.1), -0.1, 1e-5 );
  
  // Point on polygon OK
  BOOST_CHECK_SMALL(distance(1.4, 1.6), (double)0.0001);

  // Point on polygon OK
  BOOST_CHECK_SMALL(distance(1.4, 0.4), (double)0.0001);
}


BOOST_AUTO_TEST_CASE(PolygonConcaveTest)
{
  double poly[10] = {0,0,
                     0,2,
                     1,1,
                     2,2,
                     2,0};

  auto distance = [poly](const double x, const double y)
  {
    double sdf = std::numeric_limits<double>::signaling_NaN();
    getSDF(&sdf,x-0.5,y-0.5,1.0,1.0,0,0,1,1,poly,5,100);
    return sdf;
  };

  // Point (0.5,0.5) inside polygon OK
  BOOST_CHECK_GT(distance(0.5,0.5), 0);

  BOOST_CHECK_CLOSE(distance(0.5,0.5), 0.5,1e-5);

  // Point (1.5, 0.5) inside polygon
  BOOST_CHECK_GT(distance(1.5,0.5), 0);

  // Closest distance goes to vertical side
  BOOST_CHECK_CLOSE(distance(1.5,0.5), 0.5, 1e-5);
  
  // Point (0.4,1.4) inside polygon
  BOOST_CHECK_GT(distance(0.4,1.4), 0);

  // Shortest distance goes along perp. to sloping side
  BOOST_CHECK_CLOSE(distance(0.4,1.4), sqrt(2*0.1*0.1), 1e-5);
  
  // Point (1.6,1.4) inside polygon
  BOOST_CHECK_GT(distance(1.6,1.4), 0);

  // Shortest distance goes along perp. to sloping side
  BOOST_CHECK_CLOSE(distance(1.6,1.4), sqrt(2*0.1*0.1), 1e-5);
  
  // Point (1,1.5) outside polygon
  BOOST_CHECK_LT(distance(1,1.5), 0);

  // Shortest dist. goes along perp. to sloping side
  BOOST_CHECK_CLOSE(distance(1,1.5), -0.5*sqrt(2*0.5*0.5), 1e-5);

  // Point (3,1) outside polygon
  BOOST_CHECK_LT(distance(3,1), 0);

  BOOST_CHECK_CLOSE(distance(3,1), -1, 1e-5);

  // Check (0.4,1.6) on polygon
  BOOST_CHECK_SMALL(distance(0.4,1.6), (double)0.0001);

  // Check (1.6,1.6) on polygon
  BOOST_CHECK_SMALL(distance(1.6,1.6), (double)0.0001);

  // Check that shortest distance from (1,3) to polygon goes to one of the top vertices
  BOOST_CHECK_CLOSE(distance(1,3), -sqrt(2), 1e-5);

  // Check that shortest dist. from (-0.3,-0.1) goes to bottom left vertex
  BOOST_CHECK_CLOSE(distance(-0.3,-0.1), -sqrt(0.3*0.3 + 0.1*0.1), 1e-5);

  // Check that nearest point to (1,0.9) is concave point (not lower line, for example)
  BOOST_CHECK_CLOSE(distance(1,0.9), 0.1, 1e-5);
  
  // Check that dist. to pt directly in line with an edge-line is OK.
  BOOST_CHECK_CLOSE(distance(2,-0.2), -0.2, 1e-5);
}

BOOST_AUTO_TEST_SUITE_END()

