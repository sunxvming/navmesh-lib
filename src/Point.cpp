#include "Point.h"

Point::Point(double x, double y)
{
	this->x = x;
	this->y = y;
}


int32_t Point::visible(Point p, bool equal = true)	
{
	if (equal) return p.y*x <= p.x*y;
	return p.y*x < p.x*y;
}

Point::~Point()
{
}
