#include "gtest/gtest.h"
#include <iostream>
#include<ctime>
#include<vector>
#include<math.h>
#include <vector>
#include<opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/core/core.hpp>
#include <cstdio>
#include <fstream>
using namespace std;
#define N_CELLS 300
#define M_CELLS 300
#define pi 3.1415926

vector<cv::Point3d> poly;
cv::Point3d a = cv::Point3d(100.0,100.0,0);
cv::Point3d b = cv::Point3d(0,0,0);
cv::Point3d c = cv::Point3d(50.0,0,0);
extern long double DistancePointAPolygon(vector<cv::Point3d>& polyline, cv::Point3d& pt);

namespace {
TEST(DistancePointAPolygonTest, Eq) {
  EXPECT_EQ(0.0,DistancePointAPolygon(poly,a));
  EXPECT_EQ(100.0,DistancePointAPolygon(poly,b));
  EXPECT_EQ(50.0,DistancePointAPolygon(poly,c));
}
}

int main(int argc, char **argv)
{
 testing::InitGoogleTest(&argc, argv);
 poly.push_back(cv::Point3d(100.0, 100.0, 0));
 poly.push_back(cv::Point3d(-100.0, 100.0, 0));
 poly.push_back(cv::Point3d(100.0, -100.0, 0));
 poly.push_back(cv::Point3d(-100.0, -100.0, 0));
 return RUN_ALL_TESTS();
}
