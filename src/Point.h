#pragma once
#include "stdint.h"
class Point
{
public:
	double x;
	double y;
public:
	Point(double x,double y);
	Point() { this->x = 0; this -> y = 0; }
	Point(const Point& p) { this->x = p.x; this->y = p.y; }
	double Dot(Point p){ return this->x*p.x+this->y*p.y; };
	int32_t visible(Point p, bool equal );		//p是否在右侧
	bool operator == (const Point p) { return x == p.x&&y == p.y; };//这个地方就是用==，因为没有精度丢失
	Point operator - (const Point p) { return Point(this->x-p.x, this->y-p.y); };
	~Point();
};

