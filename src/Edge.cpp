#include "Edge.h"
#include "Polygon.h"
#include "log.h"

Edge::Edge(int32_t t1, int32_t t2, int32_t p1, int32_t p2)
{
	triangles[0] = t1;
	triangles[1] = t2;
	points[0] = p1;
	points[1] = p2;
}


int32_t Edge::IsRestrain(Polygon* p)
{
	int32_t p1 = points[0];
	int32_t p2 = points[1];
	int32_t dp = p2 - p1;
	if (dp == 1 || dp == -1) return 1;
	if (p1 && p2) return 0;
	return p1 + p2 == p->points.size() - 1;
}

int Edge::visible(Polygon* p, Point pt, bool equal)
{
	Point p1 = p->GetPoint(points[0]);
	Point p2 = p->GetPoint(points[1]);
	return Math::visible(p1, p2, pt, equal);
}

double Edge::Distance(Polygon* p, Point pt)
{
	Point p1 = p->GetPoint(points[0]);
	Point p2 = p->GetPoint(points[1]);
	return Math::distancePtSeg(pt, p1, p2);
}

Point Edge::FindPoint(Polygon*p, Point pt)
{
	Point p1 = p->GetPoint(points[0]);
	Point p2 = p->GetPoint(points[1]);
	
	double x2 = p2.x;
	double y2 = p2.y;
	double x1 = p1.x;
	double y1 = p1.y;
	double y3 = pt.y;
	double x3 = pt.x;

	double p12x = p2.x - p1.x;
	double p12y = p2.y - p1.y;
	double d1x = pt.x - p1.x;
	double d1y = pt.y - p1.y;
	Point inter;
	if (p12x*d1x + p12y*d1y <= 0)
	{
		inter = p1;
	}
	else {
		double d2x = pt.x - p2.x;
		double d2y = pt.y - p2.y;
		if (p12x*d2x + p12y*d2y >= 0) {
			inter = p2;
		}
		else {
			double a = x2 - x1;
			double b = y2 - y1;
			double c = -b*y3 - a*x3;
			double c2 = b*x1 + a*y1;

			double y4 = (a*c2 - b*c) / (a*a + b*b);
			double x4 = -(a*c+b*c2)/ (a*a + b*b);

			double x34 = x4 - x3;
			double y34 = y4 - y3;
			double d34 = sqrt(x34*x34 + y34*y34);
			inter = Point(x4, y4);
		}
	}
	double x34 = inter.x - x3;
	double y34 = inter.y - y3;
	double d34 = sqrt(x34*x34 + y34*y34);
	return Point(inter.x + x34 * EdgeWidth *EdgeWidth / d34, inter.y + y34 * EdgeWidth *EdgeWidth / d34);
}

Edge::~Edge()
{
}
