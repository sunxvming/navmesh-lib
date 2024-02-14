#pragma once
#include <stdint.h>
#include "Point.h"
class Polygon;
class Edge
{
public:
#if WIN32 || ANDROID || LINUX
	int32_t nothing;
#endif

	int32_t triangles[2];
	int32_t points[2];
	Edge() { this->triangles[0] = -1; this->triangles[1] = -1; this->points[0] = -1; this->points[1] = -1; }
	Edge(int32_t t1, int32_t t2, int32_t p1, int32_t p2);
	int32_t IsRestrain(Polygon* p);
	int visible(Polygon* p, Point pt, bool equal);
	double Distance(Polygon* p, Point pt);
	Point FindPoint(Polygon*p, Point pt);
	virtual ~Edge();
};

