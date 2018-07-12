#define BOOST_TEST_MODULE Main
#define BOOST_TEST_DYN_LINK

#include <cmath>
#include <cassert>
#include <iostream>
#include <limits>
#include <fstream>
#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <boost/range/iterator_range.hpp>
#include "SignedDistance.hpp"

std::vector<std::array<double,2> > getNaca0012Polygon();
std::vector<std::array<double,2> > getNacaFilmPolygon();

BOOST_AUTO_TEST_SUITE(SignedDistanceTest)

double distance1(const double x, const double y)
{
  const double poly[10] = {0,0,
                     0,2,
                     1,2,
                     2,1,
                     1,0};

  std::shared_ptr<double> sdf = std::make_shared<double>();
  getSDF(sdf.get(),x-0.5,y-0.5,1.0,1.0,1,1,poly,5,100);
  return -(*sdf);
}

BOOST_AUTO_TEST_CASE(PolygonTest)
{
  // Point inside polygon OK
  BOOST_CHECK_GT(distance1(0.5,1), 0);

  BOOST_CHECK_CLOSE(distance1(0.5,1), 0.5, 1e-5);

  // Point inside polygon OK
  BOOST_CHECK_GT(distance1(1.5,1), 0);

  // Shortest distance goes along line perp. to sloping side
  BOOST_CHECK_CLOSE(distance1(1.5,1), 0.5 * sqrt(2*0.5*0.5), 1e-5);
  
  // Point outside polygon OK
  BOOST_CHECK_LT(distance1(2,0), 0 );

  // Shortest distance goes along normal to sloping side
  BOOST_CHECK_CLOSE(distance1(2,0), - 0.5 * sqrt(2*1*1), 1e-5 );
  
  // Point outside polygon OK
  BOOST_CHECK_LT(distance1(0.5, 2.1), 0 );

  BOOST_CHECK_CLOSE(distance1(0.5, 2.1), -0.1, 1e-5 );
  
  // Point on polygon OK
  BOOST_CHECK_SMALL(distance1(1.4, 1.6), (double)0.0001);

  // Point on polygon OK
  BOOST_CHECK_SMALL(distance1(1.4, 0.4), (double)0.0001);
}


double distance2(const double x, const double y)
{
  double poly[10] = {0,0,
                     0,2,
                     1,1,
                     2,2,
                     2,0};

  std::shared_ptr<double> sdf = std::make_shared<double>();
  getSDF(sdf.get(),x-0.5,y-0.5,1.0,1.0,1,1,poly,5,100);
  return -(*sdf);
}

BOOST_AUTO_TEST_CASE(PolygonConcaveTest)
{
  // Point (0.5,0.5) inside polygon OK
  BOOST_CHECK_GT(distance2(0.5,0.5), 0);

  BOOST_CHECK_CLOSE(distance2(0.5,0.5), 0.5,1e-5);

  // Point (1.5, 0.5) inside polygon
  BOOST_CHECK_GT(distance2(1.5,0.5), 0);

  // Closest distance goes to vertical side
  BOOST_CHECK_CLOSE(distance2(1.5,0.5), 0.5, 1e-5);
  
  // Point (0.4,1.4) inside polygon
  BOOST_CHECK_GT(distance2(0.4,1.4), 0);

  // Shortest distance goes along perp. to sloping side
  BOOST_CHECK_CLOSE(distance2(0.4,1.4), sqrt(2*0.1*0.1), 1e-5);
  
  // Point (1.6,1.4) inside polygon
  BOOST_CHECK_GT(distance2(1.6,1.4), 0);

  // Shortest distance goes along perp. to sloping side
  BOOST_CHECK_CLOSE(distance2(1.6,1.4), sqrt(2*0.1*0.1), 1e-5);
  
  // Point (1,1.5) outside polygon
  BOOST_CHECK_LT(distance2(1,1.5), 0);

  // Shortest dist. goes along perp. to sloping side
  BOOST_CHECK_CLOSE(distance2(1,1.5), -0.5*sqrt(2*0.5*0.5), 1e-5);

  // Point (3,1) outside polygon
  BOOST_CHECK_LT(distance2(3,1), 0);

  BOOST_CHECK_CLOSE(distance2(3,1), -1, 1e-5);

  // Check (0.4,1.6) on polygon
  BOOST_CHECK_SMALL(distance2(0.4,1.6), (double)0.0001);

  // Check (1.6,1.6) on polygon
  BOOST_CHECK_SMALL(distance2(1.6,1.6), (double)0.0001);

  // Check that shortest distance from (1,3) to polygon goes to one of the top vertices
  BOOST_CHECK_CLOSE(distance2(1,3), -sqrt(2), 1e-5);

  // Check that shortest dist. from (-0.3,-0.1) goes to bottom left vertex
  BOOST_CHECK_CLOSE(distance2(-0.3,-0.1), -sqrt(0.3*0.3 + 0.1*0.1), 1e-5);

  // Check that nearest point to (1,0.9) is concave point (not lower line, for example)
  BOOST_CHECK_CLOSE(distance2(1,0.9), 0.1, 1e-5);
  
  // Check that dist. to pt directly in line with an edge-line is OK.
  BOOST_CHECK_CLOSE(distance2(2,-0.2), -0.2, 1e-5);
}

void test_against_data_file(const std::string& filename)
{
  const double nan = std::numeric_limits<double>::signaling_NaN();
  double xMin = nan;
  double yMin = nan;
  double dx = nan;
  double dy = nan;
  double maxDist = nan;
  int width, height;
  std::vector<double> sdf_file;
  std::vector<std::array<double,2> > vertices;

  std::ifstream file;
  file.open(filename);
  BOOST_REQUIRE(file);
  {
    std::string key, c;
    file >> c >> key >> xMin;
    file >> c >> key >> yMin;
    file >> c >> key >> dx;
    file >> c >> key >> dy;
    file >> c >> key >> width;
    file >> c >> key >> height;
    file >> c >> key >> maxDist;

    std::string item;
    file >> c >> item; // # params
    file >> item;
    assert(item.compare("Polygon "));

    file >> item; // { // Close this again so editor doesn't get confused} 
    std::string sep2;
    do
    {
      std::string sep1, parL, parR;
      std::array<double,2> v;
      file >> parL >> v[0] >> sep1 >> v[1] >> parR >> sep2;
      vertices.push_back(v);
    }
    while(sep2[0] == ',');
    file >> c >> item; //# data: 
    assert(item.compare("data: "));

    // Body
    sdf_file.resize(width*height);
    if(sdf_file.size() > 0)
    {
      while(!file.eof())
      {
        int i,j;
        double sdf;
        file >> i >> j >> sdf;
        sdf_file[i+width*j] = sdf;
      }
    }
  }

  std::vector<double> sdf_comp(width*height);

  std::ostringstream chkmsg;
  chkmsg << "Testing file: " << filename << std::endl
         << "getSDF(..., " << xMin << ", " << yMin << ", " << dx << ", " << dy << ", "
         << width << ", " << height << ", ..., " << vertices.size() << ", " << maxDist << ")";
  BOOST_TEST_CHECKPOINT(chkmsg.str());
  BOOST_TEST_MESSAGE(chkmsg.str());
  getSDF(&(sdf_comp.front()),xMin,yMin,dx,dy,width,height,&(vertices.front()[0]),vertices.size(),maxDist);
  getSDF(&(sdf_comp.front()),xMin,yMin,dx,dy,width,height,&(vertices.front()[0]),vertices.size(),maxDist);

  for(int i=0; i<width*height; i++)
  {
    BOOST_CHECK_CLOSE(sdf_comp[i]/dx, -sdf_file[i]/dx, 1e-6); // - because alo uses different sign convention
  }

}

BOOST_AUTO_TEST_CASE(PolygonFileTest)
{
  boost::filesystem::path test_data("test_data");
  BOOST_REQUIRE(boost::filesystem::is_directory(test_data));
  for(auto file : boost::make_iterator_range(
        boost::filesystem::directory_iterator(test_data), {}))
  {
    test_against_data_file(file.path().string());
  }
}

BOOST_AUTO_TEST_SUITE_END()

