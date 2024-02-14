#pragma once
#include "Point.h"
#include <math.h>
class Circle
{
private:
	Point center;
	double r;
public:
	Circle(Point p1, Point p2, Point p3);
	Point GetCenter();
	double GetR();
	virtual ~Circle();
};

