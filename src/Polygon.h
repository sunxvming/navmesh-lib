#pragma once
#include<vector>
#include"Point.h"
#include"Triangle.h"
#include"Edge.h"
#include"Circle.h"
#include"Math.h"
#include<unordered_map>
#include <stdint.h>
#include <stdio.h>

using namespace std;

#define  MAXPOINT 100000
#if _WINDOWS
	typedef unordered_map<int32_t, int32_t> Hash;
#else
	//typedef std::tr1::unordered_map<int32_t, int32_t> Hash;
	typedef unordered_map<int32_t, int32_t> Hash;
#endif
#ifndef _ASSERT
	#define _ASSERT(expr) ((void)0)
#endif

#define PIndex(p1, p2) (p1>p2?(p2*MAXPOINT+p1):(p1*MAXPOINT+p2))
#define SIndex(p1, p2, index) (p1=index/MAXPOINT, p2=index%MAXPOINT)
#define contain(p) ( lt == p || rt == p || rb == p || lb == p )

class Cell {
public:
	vector<int32_t> points;
	vector<int32_t> edges;
};

class Grid {
public:
	vector<Cell> cells;
	int32_t gride;
	double minx;
	double miny;
	double maxx;
	double maxy;
	int32_t xnum;
	int32_t ynum;
};

class Line {
public:
	Point p1;
	Point p2;
	float color[3];
};

class DijCache {
public:
	int16_t preindex; //前驱节点
	float dis; //距离
};
class Polygon
{
public:
	vector<Point> points;//外围顶点,顺时针
	vector<Triangle> triangles;
	vector<Edge> edges;
	vector<vector<DijCache>> dij;
	vector<vector<DijCache>> center;
	Hash tmpedges;
	Hash tmpcells;
	vector<int32_t> pointsnum;
	Grid grid;
private:
	inline void Gen();
	inline void CreateTriangle(Hash* eindex, int32_t p1, int32_t p2, int32_t p3);
	inline int32_t CreateEdge(Hash* eindexs, int32_t triangle, int32_t p1, int32_t p2);
	inline int32_t FindDT(Grid* grid, int32_t p1, int32_t p2 );
	inline double Distance(Point from, int32_t edge, int32_t triangle, Point to);
	inline double Distance(Point from, Point center, Point to);
	void Delaunay();
public:
	Point GetPoint(int32_t p);
	bool IsIntersect(int32_t edgepos, int32_t pa1, int32_t p1);
	bool JudgeIsVisible(int32_t pa1, int32_t p1);

	static Polygon CreateFromShort(char* cont, int32_t*len);
	//static Polygon CreateFromShort2(char* cont, int32_t*len);
	Polygon(double* pos, int32_t size);
	Polygon();		//--Polygon(char *pos)
	
	vector<Line> GetLines();
	vector<Line> GetGrideLines();
	vector<Point> GetCenters();
	int32_t Contain(double x, double y);
	int32_t IsFrist(int32_t p) { return p == 4; }
	int32_t FindTriangle(Point p);
	int32_t CheckPath(double startx, double starty, double endx, double endy);
	int32_t GetNextEposId(int32_t eposId);
	int16_t GetPreIndex(int16_t from, int16_t end, bool iscenter);
	float GetDijDis(Point from, Point to, int16_t pointfrom, int16_t pointto);
	vector<Point> FindPath(Point from, Point to, bool isturn);	//true时为中点寻路，false为拐点寻路
	Point FindNear(double x, double y);
	void GenExtData();
	void InitDijCache(FILE *file);
	void Save(FILE* file);
	virtual ~Polygon();
};

