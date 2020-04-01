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
//#define k 20//随机迭代次数

//long T;
vector<cv::Point3d> mCircle;
/*class RandomNumber {
public:
	RandomNumber() {
		srand(time(NULL) + 30 * (time(NULL) - T));
		cout << time(NULL) + 30 * (time(NULL) - T) << endl;
	}
	float get(float end) {
		return ((float)end * rand() / RAND_MAX);//end为最大值，其随机域为0~end
	}
};*/

/*cv::Point3f RandomFindPIA() {
	RandomNumber dx, dy,sita;
	float tmpx, tmpy, tmpsita,x,y,qua;
	tmpsita = sita.get(2);
	if (tmpsita <= 0.5)
		qua = 1;
	else if (tmpsita <= 1)
		qua = 2;
	else if (tmpsita <= 1.5)
		qua = 3;
	else
		qua = 4;
	tmpsita *= pi;
	tmpx = (dx.get(1)) * (1 - abs(cosf(tmpsita)));
	tmpy = (dy.get(1)) * (1 - abs(sinf(tmpsita)));
	if (qua = 1) {
		x = cosf(tmpsita) + tmpx;
		y = sinf(tmpsita) + tmpy;
	}
	else if (qua = 2) {
		x = cosf(tmpsita) - tmpx;
		y = sinf(tmpsita) + tmpy;
	}
	else if (qua = 3) {
		x = cosf(tmpsita) - tmpx;
		y = sinf(tmpsita) - tmpy;
	}
	else {
		x = cosf(tmpsita) + tmpx;
		y = sinf(tmpsita) - tmpy;
	}
	cv::Point3f pia(x, y, qua);
	return pia;
}*/
long double DistancePointAPolygon(vector<cv::Point3d>& polyline, cv::Point3d& pt);
//bool IsPointInPolygon(vector<cv::Point3f>& polyline, cv::Point3f& pt);
//int get_line_intersection(cv::Point3f& p0, cv::Point3f& p1, cv::Point3f& p2, cv::Point3f& p3);
cv::Point3d GeometryFindPIA(vector<cv::Point3d> polygon, long double bounds[]);
int  FindInscribedCircleCenter(vector<cv::Point3d> polygon, cv::Point3d& output)//寻找内切圆圆心
{
	long double bounds[4] = { 100.0,-100.0,100.0,-100.0 };//先约束正方形的范围
	/*for (int i = 0; i < polygon.size(); i++)//根据所剩的多边形构成一个新的小的正方形
	{
		cv::Point3f pt = polygon[i];
		if (pt.x < bounds[0]) bounds[0] = pt.x;//正方形左边
		if (pt.x > bounds[1]) bounds[1] = pt.x;//正方形右边
		if (pt.y < bounds[2]) bounds[2] = pt.y;//正方形上边
		if (pt.y > bounds[3]) bounds[3] = pt.y;//正方形下边
	}
	*/

	cv::Point3d point_pia;
	long double flt_tmp = FLT_MAX;
	int count = 1;
	cv::Point3d point_tmp;// = new cv::Point3f;
	while(count++)
	{
		point_tmp = GeometryFindPIA(polygon, bounds);//序列法中寻找一个候选点
		//point_tmp = RandomFindPIA();//随机法寻找一个侯选点
		// update current PIA
		point_pia.x = point_tmp.x;
		point_pia.y = point_tmp.y;

		// update the bounds
		flt_tmp = (bounds[1] - bounds[0]) / ((long double) (sqrtf(2)) * 2);
		bounds[0] = (point_pia.x - flt_tmp>-100.0)? (long double) (point_pia.x - flt_tmp):-100.0;
		bounds[1] = (point_pia.x + flt_tmp < 100.0) ? (long double) (point_pia.x + flt_tmp) :100.0;
		flt_tmp = (bounds[3] - bounds[2]) / ((long double) (sqrtf(2)) * 2);
		bounds[2] = (point_pia.y - flt_tmp > -100.0) ? (long double) (point_pia.y - flt_tmp) : -100.0;
		bounds[3] = (point_pia.y + flt_tmp < 100.0) ? (long double) (point_pia.y + flt_tmp) : 100.0;


		if (bounds[1] - bounds[0] < 0.1 || bounds[3]- bounds[2] < 0.1) break;

		//	printf("Candidate PIA: (%f,%f)\n", point_pia.x, point_pia.y);
	}


	long double tmp_distance = DistancePointAPolygon(polygon, point_pia);
	output.x= point_pia.x;
	output.y= point_pia.y;
	//output[2] = polygon[0].z;
	output.z= float(tmp_distance);
	return 1;
}
cv::Point3d GeometryFindPIA(vector<cv::Point3d> polygon, long double bounds[])//找到内切圆函数
{
	cv::Point3d pia;

	pia.x = (bounds[0] + bounds[1]) / 2;
	pia.y = (bounds[2] + bounds[3]) / 2;
	cv::Point3d tmp;


	long double increment_x = double(bounds[1] - bounds[0]) / N_CELLS;
	long double increment_y = double(bounds[3] - bounds[2]) / M_CELLS;


	long double max_distance = 0;


	int i, j;
	long double tmp_distance = -1;
	for (i = 0; i <= N_CELLS; i++)
	{

		tmp.x = bounds[0] + i * increment_x;

		for (j = 0; j <= M_CELLS; j++)
		{
			tmp.y = bounds[2] + j * increment_y;


			// compare with candidate PIA if point is in polygon
			//if (IsPointInPolygon(polygon, tmp))
			//{
			tmp_distance = DistancePointAPolygon(polygon, tmp);
			if (tmp_distance > max_distance)
			{
				max_distance = tmp_distance;
				pia.x = tmp.x;
				pia.y = tmp.y;
			}
			//}
		}
	}



	return pia;
}
/*bool IsPointInPolygon(vector<cv::Point3f>& polyline, cv::Point3f& pt)//判断点是否在多边形内
{
	int count = 0;
	cv::Point3f mark1(1000, 1000, 0);
	for (int i = 0; i < polyline.size(); i++)
	{
		int id = (i + 1) % polyline.size();
		if (get_line_intersection(polyline[i], polyline[id], pt, mark1) == 1)
		{
			count++;
		}

	}
	if (count % 2 == 1)
		return true;
	else
		return false;

}*/
/*int get_line_intersection(cv::Point3f& p0, cv::Point3f& p1, cv::Point3f& p2, cv::Point3f& p3)//判断两线段有没有交点
{

	float  s_numer, t_numer, denom, t;
	cv::Point3f s10, s32;
	s10 = p1 - p0;
	s32 = p3 - p2;


	denom = s10.x * s32.y - s32.x * s10.y;
	if (denom == 0)//平行或共线
		return 0; // Collinear
	bool denomPositive = denom > 0;

	cv::Point3f s02;
	s02 = p0 - p2;

	s_numer = s10.x * s02.y - s10.y * s02.x;
	if ((s_numer < 0) == denomPositive)//参数是大于等于0且小于等于1的，分子分母必须同号且分子小于等于分母
		return 0; // No collision


	t_numer = s32.x * s02.y - s32.y * s02.x;
	if ((t_numer < 0) == denomPositive)
		return 0; // No collision

	if (fabs(s_numer) > fabs(denom) || fabs(t_numer) > fabs(denom))
		return 0; // No collision


	return 1;
}*/
long double DistancePointAPolygon(vector<cv::Point3d>& polyline, cv::Point3d& pt)//求点到多边形的距离
{
	///1 先判断垂足位置
	//cv::Point3f mark;
	long double val, valx, valy;
	long double max_dist = -1.00000;
	valx = min(100 - pt.x, pt.x + 100);
	valy = min(100 - pt.y, pt.y + 100);
	val = min(valx, valy);
	max_dist = val;
	for (int i = 4; i < polyline.size(); i++) {
		valx = ((long double) (polyline[i].x) - pt.x) * ((long double) (polyline[i].x) - pt.x);
		valy = ((long double) (polyline[i].y) - pt.y) * ((long double) (polyline[i].y)- pt.y);
		long double dist=sqrtf(valx + valy)- (long double) (polyline[i].z);
		if (dist < max_dist)
			max_dist = dist;
	}
	/*
	for (int i = 0; i < polyline.size(); i++)
	{
		int id = (i + 1) % polyline.size();

		cv::Point3f diff = polyline[id] - polyline[i];  ///  -b a

		mark.x = -diff.y * 10000 + pt.x;
		mark.y = diff.x * 10000 + pt.y;
		double dist = 0;

		if (get_line_intersection(polyline[i], polyline[id], pt, mark) == 1)  // 有交点
		{
			cv::Point3f pt_edge = pt - polyline[i];
			val = pt_edge.x * diff.x + pt_edge.y * diff.y;
			val1 = pt_edge.x * pt_edge.x + pt_edge.y * pt_edge.y;
			if (val1 < 0.00001)
				dist = 0;
			else
			{
				val1 = sqrt(val1);
				val2 = diff.x * diff.x + diff.y * diff.y;
				val2 = sqrt(val2);
				double temp = val / (val1 * val2);
				if (temp > 1) temp = 1.0;
				if (temp < -1) temp = -1.0;

				dist = val1 * sin(acos(temp));
				// double tmep_ = val / (val2*val1);
				 // dist = sqrt(1 - temp*temp)*val1;

			}


		}
		else
		{
			cv::Point3f dist1 = pt - polyline[i];
			cv::Point3f dist2 = pt - polyline[id];
			val1 = dist1.x * dist1.x + dist1.y * dist1.y;
			val2 = dist2.x * dist2.x + dist2.y * dist2.y;
			if (val1 > val2)
			{
				val2 = sqrt(val2);
				dist = val2;
			}
			else
			{
				val1 = sqrt(val1);
				dist = val1;
			}

		}
		if (dist < max_dist)
			max_dist = dist;


	}*/
	return max_dist;
}

void mCircleSolu(int M) {
	/*vector<cv::Point3f>squra;
	squra.push_back(cv::Point3f(-1, 1, 0));
	squra.push_back(cv::Point3f(-1,-1, 0));
	squra.push_back(cv::Point3f(1, 1, 0));
	squra.push_back(cv::Point3f(1,-1, 0));*/
	long double sum = 0;
	vector<cv::Point3d> c;
	cv::Point3d circle_temp;
	int i, j;
	//RandomNumber x, y;
	for (i = 0; i < M+4; ++i)
	{
		FindInscribedCircleCenter(mCircle, circle_temp);
		mCircle.push_back(circle_temp);
		//c.push_back(Circle(x.get(2.0) - 1, y.get(2.0) - 1));
	}
	for (i = 0;i < M+4;++i)
	{
		c.push_back(mCircle.back());
		mCircle.pop_back();
	}
	cv::Point3d tmp;
	for (i = 0;i < M;++i) {
		tmp = c.back();
		cout << "第" << i + 1 << "个点的坐标是: (" << tmp.x/100<<','<< tmp.y/100<<')'<<endl;
		c.pop_back();
		sum += tmp.z * tmp.z/10000;
	}
	cout <<"r的平方和为"<< sum;
}

/*int main()
{
	//double o;
	//o = pow(1, 2);
	//T = time(NULL);
	int m;
	cout << "输入圆的数量:";
	cin >> m;
	mCircle.push_back(cv::Point3d(100.0, 100.0, 0));
	mCircle.push_back(cv::Point3d(-100.0, 100.0, 0));
	mCircle.push_back(cv::Point3d(100.0, -100.0, 0));
	mCircle.push_back(cv::Point3d(-100.0, -100.0, 0));
	mCircleSolu(m);
	return 0;
}*/
