#include "Math.h"

//求两条线段的交点
static int32_t sgn(double x)
{
	return x < -(1e-6) ? -1 : (x > 1e-6);
}

static double Cross(Point p1, Point p2, Point p3, Point p4)
{
	return (p2.x - p1.x)*(p4.y - p3.y) - (p2.y - p1.y)*(p4.x - p3.x);
}

static double Area(Point p1, Point p2, Point p3)
{
	return Cross(p1, p2, p1, p3);
}

static double fArea(Point p1, Point p2, Point p3)
{
	return fabs(Area(p1, p2, p3));
}

/*#define max( a, b ) ((a)>(b)?(a):(b))
#define min( a, b ) ((a)>(b)?(b):(a))
bool Math::Meet(Point p1, Point p2, Point p3, Point p4)
{
	return max(min(p1.x, p2.x), min(p3.x, p4.x)) <= min(max(p1.x, p2.x), max(p3.x, p4.x))
		&& max(min(p1.y, p2.y), min(p3.y, p4.y)) <= min(max(p1.y, p2.y), max(p3.y, p4.y))
		&& sgn(Cross(p3, p2, p3, p4) * Cross(p3, p4, p3, p1)) >= 0
		&& sgn(Cross(p1, p4, p1, p2) * Cross(p1, p2, p1, p3)) >= 0;
}*/


//判断直线AB是否与线段CD相交  
static bool lineIntersectSide(Point A, Point B, Point C, Point D)
{
	// A(x1, y1), B(x2, y2)的直线方程为：  
	// f(x, y) =  (y - y1) * (x1 - x2) - (x - x1) * (y1 - y2) = 0  

	double fC = (C.y - A.y) * (A.x - B.x) - (C.x - A.x) * (A.y - B.y);
	double fD = (D.y - A.y) * (A.x - B.x) - (D.x - A.x) * (A.y - B.y);
	return fC * fD <= 0;
}

#define max( a, b ) ((a)>(b)?(a):(b))
#define min( a, b ) ((a)>(b)?(b):(a))
bool Math::Meet(Point a, Point b, Point c, Point d)
{
	/*
	快速排斥：
	两个线段为对角线组成的矩形，如果这两个矩形没有重叠的部分，那么两条线段是不可能出现重叠的
	*/
	if (!(min(a.x, b.x) <= max(c.x, d.x) && min(c.y, d.y) <= max(a.y, b.y) && min(c.x, d.x) <= max(a.x, b.x) && min(a.y, b.y) <= max(c.y, d.y)))//这里的确如此，这一步是判定两矩形是否相交																																/*特别要注意一个矩形含于另一个矩形之内的情况*/
		return false;
	/*
	跨立实验：
	如果两条线段相交，那么必须跨立，就是以一条线段为标准，另一条线段的两端点一定在这条线段的两段
	也就是说a b两点在线段cd的两端，c d两点在线段ab的两端
	*/
	double u, v, w, z;//分别记录两个向量
	u = (c.x - a.x)*(b.y - a.y) - (b.x - a.x)*(c.y - a.y);
	v = (d.x - a.x)*(b.y - a.y) - (b.x - a.x)*(d.y - a.y);
	w = (a.x - c.x)*(d.y - c.y) - (d.x - c.x)*(a.y - c.y);
	z = (b.x - c.x)*(d.y - c.y) - (d.x - c.x)*(b.y - c.y);
	return (u*v <= 0 && w*z <= 0);
}

Point Math::Inter(Point p1, Point p2, Point p3, Point p4)
{
	double a12 = p1.y - p2.y;
	double b12 = p2.x - p1.x;
	double c12 = p1.x*p2.y -p2.x*p1.y;

	double a34 = p3.y - p4.y;
	double b34 = p4.x - p3.x;
	double c34 = p3.x*p4.y - p4.x*p3.y;
	double D = a12*b34 - a34*b12;
	Point p;
	p.x = (b12*c34 -b34*c12) / D;
	p.y = (c12*a34 -c34*a12) / D;
	return p;

	//double s1 = fArea(p1, p2, p3), s2 = fArea(p1, p2, p4);
	//return Point((p4.x*s1 + p3.x*s2) / (s1 + s2), (p4.y*s1 + p3.y*s2) / (s1 + s2));
}

int32_t Math::visible(Point p1, Point p2, Point p3, bool equal)			//p3是否在p1p2的右侧
{
	double x1 = p1.x, y1 = p1.y;
	double x2 = p2.x, y2 = p2.y;
	double x3 = p3.x, y3 = p3.y;
	if(equal) return (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3) <= 0;
	return (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3) < 0;
}

/*int32_t Math::visible(Point p2, Point p3, bool equal)		//p2是否在p1p3的右侧,equal，表示共线算不算/p2是否
{
	double x2 = p2.x, y2 = p2.y;
	double x3 = p3.x, y3 = p3.y;
	if (equal) return y3*x2 >= x3*y2;
	return y3*x2 > x3*y2;
}*/

double Math::distancePtSeg(const Point pt, const Point p, const Point q)
{
	double pqx = q.x - p.x;
	double pqy = q.y - p.y;
	double dx = pt.x - p.x;
	double dy = pt.y - p.y;
	double d = pqx*pqx + pqy*pqy;      // qp线段长度的平方
	double t = pqx*dx + pqy*dy;         // p pt向量 点积 pq 向量（p相当于A点，q相当于B点，pt相当于P点）
	if (d > 0)         // 除数不能为0; 如果为零 t应该也为零。下面计算结果仍然成立。                   
		t /= d;    // 此时t 相当于 上述推导中的 r。
	if (t < 0)
		t = 0;     // 当t（r）< 0时，最短距离即为 pt点 和 p点（A点和P点）之间的距离。
	else if (t > 1)
		t = 1;     // 当t（r）> 1时，最短距离即为 pt点 和 q点（B点和P点）之间的距离。
				// t = 0，计算 pt点 和 p点的距离; t = 1, 计算 pt点 和 q点 的距离; 否则计算 pt点 和 投影点 的距离。
	dx = p.x + t*pqx - pt.x;
	dy = p.y + t*pqy - pt.y;
	return dx*dx + dy*dy;
}

double Math::pointmul(Point from, Point p1, Point p2)
{
	return (p1.x - from.x) * (p2.x - p1.x) + (p1.y - from.y) * (p2.y - p1.y);
}

Point Math::pendalpoint(Point from, Point p1, Point p2) //from到p1 p2的垂点
{
	if ((p2.y - p1.y) * (from.x - p1.x) == (p2.x - p1.x) * (from.y - p1.y))
		return from;
	Point p;
	double x1 = p1.x;
	double x2 = p2.x;
	double x3 = from.x;
	double y1 = p1.y;
	double y2 = p2.y;
	double y3 = from.y;
	p.x = (((y2 - y1)* y3 + (x2 - x1) * x3) * (x2 - x1) - (x2 * y1 - x1 * y2) * (y2 - y1)) / ((y2 - y1) * (y2 - y1) +
			(x2 - x1) * (x2 - x1));
	p.y = ((x2 - x1)* (x2 * y1 - x1 * y2) + (y2 - y1) * ((y2 - y1) * y3 + (x2 - x1) * x3)) / ((y2 - y1) * (y2 - y1) +
			(x2 - x1) * (x2 - x1));
	return p;
}