#pragma once
#include"Polygon.h"
#include <stdint.h>

class Math {
public:
	static bool Meet(Point p1, Point p2, Point p3, Point p4);
	static Point Inter(Point p1, Point p2, Point p3, Point p4);
	static int32_t visible(Point p1, Point p2, Point p3, bool equal = true);			//p3是否在p1p2的右侧
	//static int32_t visible(Point p2, Point p3, bool equal = true);		//p2是否在p1p3的右侧,equal，表示共线算不算/p2是否
	static double distancePtSeg(Point pt, Point p, Point q);
	static double pointmul(Point from, Point p1, Point p2);
	static Point pendalpoint(Point from, Point p1, Point p2);
};