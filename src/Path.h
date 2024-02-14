#pragma once

#include<vector>
#include"Polygon.h"
#include"Math.h"
#include<math.h>
#include <stdint.h>


using namespace std;
class Polygon;
class Path
{
public:
	vector<Polygon> polygons;
	vector<double> points;
	vector<int32_t > indexs;
	vector<double> finalpath;
	vector<double> nearpoint;
	vector<double> cropoint;
	vector<int32_t> sizes;
	Hash edgehashs;
private:
	//inline char* readFile(const char* file, int32_t* size);
public:
	Path();
	//Path(const char* file);
	Path(char* cont, int32_t size);
	const double* GetPoints(int32_t* length);
	const int32_t* GetIndexs(int32_t* length);
	const int32_t* GetSizes(int32_t* length);
	bool IsInPolygon(Point pos);
	const double* FindPaths(Point start,Point end, bool isturn, int32_t* size);
	const double* FindCross(double startx, double starty, double facex, double facey);
	const double* FindNear(double x, double y);
	int32_t GetNextEposId(int32_t eposId, Polygon* p);
	void Save(const char* filename);
	void Load(const char* filename);
	int32_t CheckPath(double startx, double starty, double endx, double endy);
	~Path();
};

