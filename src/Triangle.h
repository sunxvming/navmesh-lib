#pragma once
#include "Point.h"
#include <stdint.h>

class Polygon;
#pragma pack(push)
class Triangle
{
private:
#if WIN32 || ANDROID || LINUX
	int32_t nothing;
#endif
public:
	int32_t p1;
	int32_t p2;
	int32_t p3;
	int32_t edges[3];
	Point icenter;//опл─
	Point lt;
	Point rb;
public:
	Triangle(int32_t p1, int32_t p2, int32_t p3);
	Triangle() { 
		this->p1 = 0; 
		this->p2 = 0; 
		this->p3 = 0; 
		this->edges[0] = 0;
		this->edges[1] = 0;
		this->edges[2] = 0;
	}
	bool InCircle(Polygon* polygon, int32_t p);
	void GenExtData(Polygon* p);
	int32_t Contain(Polygon* p, Point pt);
	virtual ~Triangle();
};
#pragma pack(4)
#pragma pack(pop)