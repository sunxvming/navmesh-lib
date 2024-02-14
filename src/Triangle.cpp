#include "Triangle.h"
#include"Polygon.h"
#include <string.h>

Triangle::Triangle(int32_t p1, int32_t p2, int32_t p3):p1(p1),p2(p2),p3(p3)
{
	memset(edges, -1, 3);
}

void  circle_center(Point  *center, Point  pt[3], double  *radiu)
{
	double  x1, x2, x3, y1, y2, y3;
	double  x = 0;
	double  y = 0;

	x1 = pt[0].x;
	x2 = pt[1].x;
	x3 = pt[2].x;
	y1 = pt[0].y;
	y2 = pt[1].y;
	y3 = pt[2].y;

	x = ((y2 - y1)*(y3*y3 - y1*y1 + x3*x3 - x1*x1) - (y3 - y1)*(y2*y2 - y1*y1 + x2*x2 - x1*x1)) / (2 * (x3 - x1)*(y2 - y1) - 2 * ((x2 - x1)*(y3 - y1)));
	y = ((x2 - x1)*(x3*x3 - x1*x1 + y3*y3 - y1*y1) - (x3 - x1)*(x2*x2 - x1*x1 + y2*y2 - y1*y1)) / (2 * (y3 - y1)*(x2 - x1) - 2 * ((y2 - y1)*(x3 - x1)));

	center->x = x;
	center->y = y;
	*radiu = (pt[0].x - x)*(pt[0].x - x) + (pt[0].y - y)*(pt[0].y - y);
}

// 求 aij 的代数余子式
vector< vector<double> > Cofactor(vector< vector<double> > vecDet_ij, int32_t i, int32_t j)
{
	int32_t k;
	vector< vector<double> > vecReturn;
	vector< vector<double> >::iterator veck;
	// vector<double>::iterator vecl;

	// 初始化二维容器 vecReturn
	k = 0;
	for (veck = vecDet_ij.begin(); veck<vecDet_ij.end(); veck++)
	{
		if ((veck - vecDet_ij.begin()) != i)
		{
			vecReturn.push_back(*veck); // 加入除第 i 行外的所有行

			vecReturn[k].erase(vecReturn[k].begin() + j);
			k++;
		}
	}
	return vecReturn;
}
// 计算行列式的值, 采用递归
double det_Array(vector< vector<double> > vecDet)
{
	int32_t i;
	double Sum = 0.0;
	vector< vector<double> >::iterator vec_Row = vecDet.begin();
	vector<double>::iterator vec_Column = (*vec_Row).begin();
	if (vecDet.size() == 1)
	{
		return *vec_Column;
	}
	else
	{
		i = 0;
		for (; vec_Column < (*vec_Row).end(); vec_Column++)
		{
			Sum = Sum + pow(-1.0, i) * (*vec_Column) * det_Array(Cofactor(vecDet, 0, i++));
		}
		return Sum;
	}
}

bool Triangle :: InCircle(Polygon* polygon, int32_t p)
{
	Point Vrtx0 = polygon->GetPoint(p1);
	Point Vrtx1 = polygon->GetPoint(p2);
	Point Vrtx2 = polygon->GetPoint(p3);
	Point Vrtx = polygon->GetPoint(p);

	double Radius_2;  // 半径的平方
					  /*    double Radius_2T1, Radius_2T2;*/
	double Cntrx, Cntry; //圆心坐标

						 // 求圆心和半径
	vector< vector<double> > vecDet1, vecDet2;
	vector<double> A, B, C;

	// 计算圆心的 x 坐标
	A.push_back(pow(Vrtx0.x, 2.0) + pow(Vrtx0.y, 2.0)); A.push_back(Vrtx0.y); A.push_back(1.0);
	B.push_back(pow(Vrtx1.x, 2.0) + pow(Vrtx1.y, 2.0)); B.push_back(Vrtx1.y); B.push_back(1.0);
	C.push_back(pow(Vrtx2.x, 2.0) + pow(Vrtx2.y, 2.0)); C.push_back(Vrtx2.y); C.push_back(1.0);
	vecDet1.push_back(A); vecDet1.push_back(B); vecDet1.push_back(C);

	A.clear(); B.clear(); C.clear();
	A.push_back(Vrtx0.x); A.push_back(Vrtx0.y); A.push_back(1.0);
	B.push_back(Vrtx1.x); B.push_back(Vrtx1.y); B.push_back(1.0);
	C.push_back(Vrtx2.x); C.push_back(Vrtx2.y); C.push_back(1.0);
	vecDet2.push_back(A); vecDet2.push_back(B); vecDet2.push_back(C);
	Cntrx = det_Array(vecDet1) / (2 * det_Array(vecDet2));

	// 计算圆心的 y 坐标
	vecDet1.clear(); vecDet2.clear();
	A.clear(); B.clear(); C.clear();
	A.push_back(Vrtx0.x); A.push_back(pow(Vrtx0.x, 2.0) + pow(Vrtx0.y, 2.0)); A.push_back(1.0);
	B.push_back(Vrtx1.x); B.push_back(pow(Vrtx1.x, 2.0) + pow(Vrtx1.y, 2.0)); B.push_back(1.0);
	C.push_back(Vrtx2.x); C.push_back(pow(Vrtx2.x, 2.0) + pow(Vrtx2.y, 2.0)); C.push_back(1.0);
	vecDet1.push_back(A); vecDet1.push_back(B); vecDet1.push_back(C);

	A.clear(); B.clear(); C.clear();
	A.push_back(Vrtx0.x); A.push_back(Vrtx0.y); A.push_back(1.0);
	B.push_back(Vrtx1.x); B.push_back(Vrtx1.y); B.push_back(1.0);
	C.push_back(Vrtx2.x); C.push_back(Vrtx2.y); C.push_back(1.0);
	vecDet2.push_back(A); vecDet2.push_back(B); vecDet2.push_back(C);
	Cntry = det_Array(vecDet1) / (2 * det_Array(vecDet2));

	// 外接圆的半径的平方
	Radius_2 = pow(Vrtx0.x - Cntrx, 2.0) + pow(Vrtx0.y - Cntry, 2.0);

	// 判断 Vrtx0 是否在外接圆内或其上，若是，返回 true,否则，返回 false
	double Rad_V0Cntr;
	Rad_V0Cntr = pow(Vrtx.x - Cntrx, 2.0) + pow(Vrtx.y - Cntry, 2.0);
	return Rad_V0Cntr <= Radius_2;
}

#define max( a, b ) ((a)>(b)?(a):(b))
#define min( a, b ) ((a)>(b)?(b):(a))
void Triangle::GenExtData(Polygon* p)
{
	Point pt1 = p->GetPoint(p1);
	Point pt2 = p->GetPoint(p2);
	Point pt3 = p->GetPoint(p3);
	
	icenter.x = (pt1.x+pt2.x+pt3.x)/3;
	icenter.y = (pt1.y+pt2.y+pt3.y)/3;

	double maxx = max(pt1.x, pt2.x);
	double maxy =  max(pt1.y, pt2.y);
	double minx =  min(pt1.x, pt2.x);
	double miny =  min(pt1.y, pt2.y);

	if(pt3.x>maxx) maxx = pt3.x;
	if(pt3.y>maxy) maxy = pt3.y;
	if(pt3.x<minx) minx = pt3.x;
	if(pt3.y<miny) miny = pt3.y;

	lt.x = minx, lt.y = miny;
	rb.x = maxx, rb.y = maxy;
}

int32_t Triangle::Contain(Polygon* p, Point pt)
{
	double x = pt.x, y = pt.y;
	if(x<lt.x) return 0;
	if(x>rb.x) return 0;
	if(y<lt.y) return 0;
	if(y>rb.y) return 0;

	Point pt1 = p->GetPoint(p1);
	Point pt2 = p->GetPoint(p2);
	Point pt3 = p->GetPoint(p3);

	Point v0 = pt3-pt1;
	Point v1 = pt2-pt1;
	Point v2 = pt-pt1;

	double dot00 = v0.Dot(v0) ;
    double dot01 = v0.Dot(v1) ;
    double dot02 = v0.Dot(v2) ;
    double dot11 = v1.Dot(v1) ;
    double dot12 = v1.Dot(v2) ; 
	double inverDeno = 1 / (dot00 * dot11 - dot01 * dot01) ;
    double u = (dot11 * dot02 - dot01 * dot12) * inverDeno ;
    if (u < 0 || u > 1) // if u out of range, return directly
    {
        return 0 ;
    }

    double v = (dot00 * dot12 - dot01 * dot02) * inverDeno ;
    if (v < 0 || v > 1) // if v out of range, return directly
    {
        return 0 ;
    }

    return u + v <= 1 ? 1:0 ;
}

Triangle::~Triangle()
{
}
