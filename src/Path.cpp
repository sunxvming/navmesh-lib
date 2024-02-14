#define _CRT_SECURE_NO_WARNINGS
#include "Path.h"
#include "log.h"
#include <errno.h>
#include<string.h>
Path::Path()
{	
}

Path::Path(char* cont, int32_t size)
{
	int32_t len = 0;
	while (size > 0)
	{
		//printf("sizenow %d %d\n", size, polygons.size());
		polygons.push_back(Polygon::CreateFromShort(cont, &len));
		cont += len;
		size -= len;
		//printf("polygosss %d\n", size);
	}
}

void Path::Save(const char* filename)
{
	FILE* file;
	if ((file = fopen(filename, "wb")) == 0)
		return;
	int32_t len = polygons.size();
	fwrite(&len, sizeof(int32_t), 1, file);
	for (int32_t i = 0; i < len; i++)
	{
		polygons[i].Save(file);
	}
	fclose(file);
}

static inline long filepos(FILE* f)
{
	fpos_t pos;
	if (!fgetpos(f, &pos))
#if LINUX
		return pos.__pos;
#else
		return pos;
#endif
	return -1;
}

void Path::Load(const char* filename) {
	FILE* file;
	LOGD("Path::Load %s step0", filename);
	if ((file = fopen(filename, "rb")) == 0)
	{
		LOGD("Path::Loadfail %s ", strerror(errno));
		return;
	}
		
	//16 24 80 4 8
	LOGD( "Path::Load %s step1 %d %d %d %d %d", filename, sizeof(Point), sizeof(Edge), sizeof(Triangle), sizeof(int32_t), sizeof(double));
	int32_t len;
	fread(&len, 1, sizeof(int32_t), file);
	LOGD("Path::Load %s step2 %d, offset %d", filename, len, filepos(file));
	for (int32_t i = 0; i < len; i++) {
		Polygon polygon;
		int32_t pointsLen;
		fread(&pointsLen, 1, sizeof(int32_t), file);
		//int32_t ** pr = new int32_t*[pointsLen];
		LOGD("Path::Load %s step3 %d, offset %d", filename, pointsLen, filepos(file));
		for (int32_t j = 0; j < pointsLen; j++) {
			Point point;
			fread(&(point), 1, sizeof(Point), file);
			polygon.points.push_back(point);
		}

		int32_t triangeLen;
		fread(&triangeLen, 1, sizeof(int32_t), file);
		LOGD("Path::Load %s step4 %d, offset %d", filename, triangeLen, filepos(file));
		for (int32_t j = 0; j < triangeLen; j++) {
			Triangle tri;
			fread(&tri, 1, sizeof(Triangle), file);
			polygon.triangles.push_back(tri);
		}

		int32_t edgeLen;
		fread(&edgeLen, 1, sizeof(int32_t), file);
		LOGD("Path::Load %s step5 %d, offset %d", filename, edgeLen, filepos(file));
		for (int32_t j = 0; j < edgeLen; j++) {
			Edge edge;
			fread(&(edge), 1, sizeof(Edge), file);
			polygon.edges.push_back(edge);
		}

		int32_t pLen;
		fread(&pLen, 1, sizeof(int32_t), file);
		LOGD("Path::Load %s step6 %d, offset %d", filename, pLen, filepos(file));
		for (int32_t j = 0; j < pLen; j++) {
			int32_t num;
			fread(&(num), 1, sizeof(int32_t), file);
			polygon.pointsnum.push_back(num);
		}
	
		int32_t cellsize;
		fread(&cellsize, 1, sizeof(int32_t), file);
		Grid grid;
		LOGD("Path::Load %s step7 %d, offset %d", filename, cellsize, filepos(file));
		for (int32_t j = 0; j < cellsize; j++)
		{
			Cell cel;
			int32_t gridpoLen;
			fread(&gridpoLen, 1, sizeof(int32_t), file);
			for (int32_t k = 0; k < gridpoLen; k++) {
				int32_t pk;
				fread(&(pk), 1, sizeof(int32_t), file);
				cel.points.push_back(pk);
			}

			int32_t gridedLen;
			fread(&gridedLen, 1, sizeof(int32_t), file);
			for (int32_t k = 0; k < gridedLen; k++) {
				int32_t ek;
				fread(&(ek), 1, sizeof(int32_t), file);
				cel.edges.push_back(ek);
			}
			grid.cells.push_back(cel);
		}
		int32_t gride;
		fread(&(gride), 1, sizeof(int32_t), file);
		grid.gride = gride;

		double minx, maxx, miny, maxy;
		fread(&(minx), 1, sizeof(double), file);
		grid.minx = minx;
		fread(&(maxx), 1, sizeof(double), file);
		grid.maxx = maxx;
		fread(&(miny), 1, sizeof(double), file);
		grid.miny = miny;
		fread(&(maxy), 1, sizeof(double), file);
		grid.maxy = maxy;
		//printf("%f %f %f %f==path max \n", minx, maxx, miny, maxy);
		int32_t xnum, ynum;
		fread(&(xnum), 1, sizeof(int32_t), file);
		grid.xnum = xnum;
		fread(&(ynum), 1, sizeof(int32_t), file);
		grid.ynum = ynum;
		
		polygon.grid = grid;
		//polygon.InitDijCache();
		int32_t size = polygon.points.size();
		//polygon.dij = new DijCache *[size];
		int16_t preindex;
		float disnow;
		polygon.dij.resize(size);
		for (int32_t i = 0; i < size; i++)
		{
			polygon.dij[i].resize(size);
		}
		for(int32_t i = 0; i < size; i++)
		{
			for (int32_t j = i; j < size; j++)
			{
				if (i == j)
				{
					polygon.dij[i][j].preindex = i;
					polygon.dij[i][j].dis = 0;
					continue;
				}
				fread(&(preindex), 1, sizeof(int16_t), file);
				//fread(&(backpreindex), 1, sizeof(int16_t), file);
				fread(&(disnow), 1, sizeof(float), file);
				//if (i == 14 && j == 243)
					//printf("polygon.dij[i][j] %d %d %f %f %d\n", i, j, disnow, polygon.dij[i][j].dis, preindex);
				polygon.dij[i][j].preindex = preindex;
				//polygon.dij[j][i].preindex = backpreindex;
				polygon.dij[i][j].dis = disnow;
				//polygon.dij[j][i].dis = disnow;
				//printf("polygon.dij[i][j] %d %d %f %f %d %d\n", i, j, disnow, polygon.dij[i][j].dis, preindex, backpreindex);
			}
		}
		size = polygon.triangles.size();
		polygon.center.resize(size);
		for (int32_t i = 0; i < size; i++)
		{
			polygon.center[i].resize(size);
		}
		for (int32_t i = 0; i < size; i++)
		{
			for (int32_t j = i; j < size; j++)
			{
				if (i == j)
				{
					polygon.center[i][j].preindex = i;
					polygon.center[i][j].dis = 0;
					continue;
				}
				fread(&(preindex), 1, sizeof(int16_t), file);
				//fread(&(backpreindex), 1, sizeof(int16_t), file);
				fread(&(disnow), 1, sizeof(float), file);
				polygon.center[i][j].preindex = preindex;
				//polygon.center[j][i].preindex = backpreindex;
				polygon.center[i][j].dis = disnow;
				//polygon.center[j][i].dis = disnow;
				//printf("polygon.dij[i][j] %d %d %f %f %d %d\n", i, j, disnow, polygon.dij[i][j].dis, preindex, backpreindex);
			}
		}
		//printf("init %d", polygons.size());
		polygons.push_back(polygon);
	}
	LOGD("Path::Load %s step8, offset %d", filename, filepos(file));
	fclose(file);
}

const int32_t* Path::GetSizes(int32_t* length) {
	sizes.clear();
	for (int32_t i = 0; i < polygons.size(); i++) {
		const vector<int32_t> & pointsnum = polygons[i].pointsnum;
		sizes.insert(sizes.end(), pointsnum.begin( ), pointsnum.end());
	}
	*length = sizes.size();
	return sizes.data();
}

const double* Path::GetPoints(int32_t* length) {
	points.clear();

	for (int32_t i = 0; i < polygons.size(); i++) {
		for (int32_t j = 0; j < polygons[i].points.size(); j++) {
			Point p = polygons[i].points[j];
			points.push_back(p.x);
			points.push_back(p.y);
		}
	}
	*length = points.size();
	return points.data();
}

const int32_t* Path::GetIndexs(int32_t* length) {
	indexs.clear();
	int32_t vector = 0;
	for (int32_t i = 0; i < polygons.size(); i++) {
		for (int32_t j = 0; j < polygons[i].triangles.size(); j++)
		{
			indexs.push_back(polygons[i].triangles[j].p1 + vector);
			indexs.push_back(polygons[i].triangles[j].p2 + vector);
			indexs.push_back(polygons[i].triangles[j].p3 + vector);
		}
		vector += polygons[i].points.size();
	}
	*length = indexs.size();
	return indexs.data();
}

const double* Path::FindPaths(Point start, Point end, bool isturn,int32_t* size) {
	//printf("Path::FindPaths============");
	finalpath.clear();
	int32_t startPIndex = -1;
	int32_t endPIndex = -1;
	//printf("Path::FindPaths=====00000000000000=======%d \n", polygons.size());
	for (int32_t i = 0; i < polygons.size(); i++) {
		
		double minx = polygons[i].grid.minx;
		double maxx = polygons[i].grid.maxx;
		double miny = polygons[i].grid.miny;
		double maxy = polygons[i].grid.maxy;
		//printf("polygons===for===%d %f %f %f %f %f %f %f %f \n", i, minx, maxx, miny, maxy, start.x, start.y, end.x, end.y);
		if ((start.x >= minx &&start.x <= maxx) && (start.y >= miny &&start.y <= maxy)) {
			startPIndex = i;
			int32_t startTIndex = polygons[startPIndex].FindTriangle(start);
		}
		if ((end.x >= minx &&end.x <= maxx) && (end.y >= miny &&end.y <= maxy)) {
			endPIndex = i;
			int32_t endTIndex = polygons[endPIndex].FindTriangle(end);
		}
		if ((startPIndex >= 0 && endPIndex >= 0 && startPIndex == endPIndex))
			break;
	}
	//printf("Path::FindPaths=====111111111111111111=======%d %d %d \n", polygons.size(), startPIndex, endPIndex);
	if ((startPIndex != endPIndex) || (startPIndex < 0) || (endPIndex < 0)) {
		//printf("Path::FindPaths===~~~~~~====%d %d ", startPIndex, endPIndex);
		*size = 0;
		//printf("Path::FindPaths====2222222222222========%d %d ", startPIndex, endPIndex);
		return 0;
	}
	else {
		int32_t pIndex = startPIndex;

		int32_t startTIndex = polygons[pIndex].FindTriangle(start);
		//printf("start triangle!!!!!!!!!! %d\n", startTIndex);
		if (startTIndex < 0) return 0;

		int32_t endTIndex = polygons[pIndex].FindTriangle(end);
		//printf("end triangle@@@@@@ %d\n", endTIndex);
		if (endTIndex < 0) return 0;

		vector<Point> f = polygons[pIndex].FindPath(start, end, isturn);
		int32_t fsize = f.size();
		for (int32_t i = fsize - 1 ; i >= 0; i--) {
			finalpath.push_back(f[i].x);
			finalpath.push_back(f[i].y);
		}
		*size = fsize * 2;
		//printf("fsize %d \n", fsize);
		if (fsize < 1)
			return 0;
		return finalpath.data();
	}
	*size = 0;
	return 0;
}



int32_t Path::GetNextEposId(int32_t eposId, Polygon* p) {
	int32_t sum = 0;
	int32_t nextId = 0;
	for (int32_t i = 0; i < p->pointsnum.size(); i++)
	{
		sum += p->pointsnum[i];
		if (eposId == sum - 1) {
			nextId = eposId - (p->pointsnum[i] - 1);
			break;
		}
		else if (eposId < sum - 1) {
			nextId = eposId + 1;
			break;
		}
		else if (eposId > sum - 1) {
			continue;
		}
	}
	return nextId;
}

bool Path::IsInPolygon(Point pos)
{
	int32_t sIndex = -1;
	int32_t stIndex = -1;
	for (int32_t i = 0; i < polygons.size(); i++) {
		double minx = polygons[i].grid.minx;
		double maxx = polygons[i].grid.maxx;
		double miny = polygons[i].grid.miny;
		double maxy = polygons[i].grid.maxy;
		if ((pos.x >= minx && pos.x <= maxx) && (pos.y >= miny &&pos.y <= maxy)) {
			sIndex = i;
			stIndex = polygons[sIndex].FindTriangle(pos);
		}
		if (sIndex >= 0 && stIndex >= 0)
			break;
	}
	if ( stIndex < 0 )
		return false;
	return true;
}

const double* Path::FindCross(double startx, double starty, double facex, double facey) {
	//startx = 135.40393315633, starty = 7.0799863646626, facex = -21.173290252686, facey = 72.436767578125;
	LOGD("Path::FindCross %f %f %f %f",  startx, starty, facex, facey);
	cropoint.clear();
	//如果facex与facey同为0，则忽略，直接返回start点
	if (facex == 0 && facey == 0) {
		cropoint.push_back(startx);
		cropoint.push_back(starty);
		return cropoint.data();
	}

	Point start = Point(startx, starty);
	Point face = Point(startx + facex, starty + facey);
	double length = 0;
	Point lpoint;
	LOGD("Path::FindCross %d, %f %f %f %f", polygons.size(), startx, starty, facex, facey);
	for (int32_t k = 0; k < polygons.size(); k++) {
		//Point vector = Point(facex - startx, facey - starty);		//向量
		double minx, maxx, miny, maxy;
		Polygon* p = &polygons[k];
		Grid* grid = &p->grid;
		minx = grid->minx;
		maxx = grid->maxx;
		miny = grid->miny;
		maxy = grid->maxy;
		LOGD("Path::FindCross1 %d  %f %f %f %f", k, minx, maxx, miny, maxy);
		if (startx<minx || startx>maxx || starty<miny || starty>maxy) continue;
		LOGD("Path::FindCross2 %d  %f %f %f %f", k, startx, starty, facex, facey);
		int32_t gride = grid->gride;
		int32_t xnum = grid->xnum;
		int32_t ynum = grid->ynum;

		//TODO::通过start face找到与边界相交的点
		Point end;
		Point lt = Point(minx, miny);
		Point lb = Point(minx, maxy);
		Point rt = Point(maxx, miny);
		Point rb = Point(maxx, maxy);

		if (facex == 0 && facey > 0) {
			Point xmax = Math::Inter(start, face, lb, rb);
			end = Point(xmax.x, xmax.y + facey);
		}
		else if (facex == 0 && facey < 0) {
			Point xmin = Math::Inter(start, face, lt, rt);
			end = Point(xmin.x, xmin.y + facey);
		}
		else if (facex > 0) {
			Point xmax = Math::Inter(start, face, rt, rb);
			end = Point(xmax.x + facex, xmax.y + facey);
		}
		else if (facex < 0) {
			Point xmin = Math::Inter(start, face, lt, lb);
			end = Point(xmin.x + facex, xmin.y + facey);
		}

		int32_t xn1 = (int32_t)((start.x - minx) / gride);
		int32_t xn2 = (int32_t)((end.x - minx) / gride);
		int32_t yn1 = (int32_t)((start.y - miny) / gride);
		int32_t yn2 = (int32_t)((end.y - miny) / gride);

		edgehashs.clear();
		double y = start.y;
		double x = start.x;
		int32_t stepx = xn2 > xn1 ? 1 : -1;
		int find = 0;
		int32_t i = xn1;
		//if (xn2 >= xn1 && xn1 < 0)i = 0;
		//if (xn1 <= xn2 && xn1 >xnum)i = xnum;
		while( 1 )
		{
			if (i > xnum || i < 0) break;
			int32_t next_y = yn2;
			double x3, y3;
			if (stepx < 0)
				x3 = i == xn2 ? end.x : i * gride + minx;
			else
				x3 = i == xn2 ? end.x : (i + stepx) * gride + minx;
			if (facex != 0)
				y3 = (x3 - start.x)*facey / facex + start.y;
			else
				y3 = starty;
			next_y = (int32_t)((y3 - miny) / gride);
			//int32_t cur_x = (int32_t)((x - minx) / gride);
			int32_t cur_y = (int32_t)((y - miny) / gride);
			int32_t stepy = next_y > cur_y ? 1 : -1;
			int32_t j = cur_y;
			while( 1 )
			{
				if (j > ynum || j < 0) break;
				int32_t edgepos = (i >= xnum ? xnum - 1 : i) + (j >= ynum ? ynum - 1 : j)*xnum;
				//if (edgepos < 0 || edgepos > xnum * ynum) continue;
				LOGD("Path::FindCross2 %d %d %d %f %f", i, j, grid->cells[edgepos].edges.size(), end.x, end.y);
				for (int32_t k = 0; k < grid->cells[edgepos].edges.size(); k++) {
					int32_t hashk = grid->cells[edgepos].edges[k];
					if (edgehashs.find(hashk) == edgehashs.end()) {
						Point p1 = p->points[hashk];
						Point p2 = p->points[GetNextEposId(hashk, p)];
						edgehashs.insert(make_pair(hashk, edgepos));
						LOGD("Path::FindCross4 %d %f %f %f %f %f %f %f %f", hashk, p1.x, p1.y, p2.x, p2.y, start.x, start.y, end.x, end.y);
						if (Math::visible(p1, p2, start, false)) continue;
						LOGD("Path::FindCross5 %d %f %f %f %f %f %f %f %f", hashk, p1.x, p1.y, p2.x, p2.y, start.x, start.y, end.x, end.y);
						if(!Math::Meet(start, end, p1, p2))  continue;
						LOGD("Path::FindCross6 %d %f %f %f %f %f %f %f %f", hashk, p1.x, p1.y, p2.x, p2.y, start.x, start.y, end.x, end.y);
						Point intep = Math::Inter(start, end, p1, p2);
						//如果是第一次进来，则首先先给lpoint和length附值
						if (length == 0) {
							lpoint = intep;
							length = (intep.x - startx)*(intep.x - startx) + (intep.y - starty)*(intep.y - starty);
						}
						//TODO:求start->face 到交点的距离 
						double curLen = (intep.x - startx)*(intep.x - startx) + (intep.y - starty)*(intep.y - starty);
						//TODO然后与length 做比较，如果更小，则记录当前的交点，否则继续下一轮。
						if (curLen < length) {
							lpoint = intep;
							length = curLen;
						}
					}
				}
				if (length > 0 && ((int32_t)((lpoint.x - minx) / gride) == i) && ((int32_t)((lpoint.y - miny) / gride) == j))
				{
					find = 1;
					break;
				}
				if (j == next_y) break;
				j += stepy;
			}
			if (find) break;
			x = x3;
			y = y3;
			if (i == xn2) break;
			i += stepx;
		}	
	}
	Point v1 = Point(lpoint.x - startx, lpoint.y - starty);
	if (length <= EdgeWidth*EdgeWidth) {
		return 0;
	}
	else {
		double l = sqrt(length);
		cropoint.push_back(lpoint.x - v1.x * 0.1 / l);
		cropoint.push_back(lpoint.y - v1.y * 0.1 / l);
	}
	return cropoint.data();
}

const double* Path::FindNear(double x, double y)
{
	nearpoint.clear();
	nearpoint.push_back(0);
	nearpoint.push_back(0);
	double dis = -1;
	for (int32_t i = 0; i < polygons.size(); i++) {
		if (!polygons[i].Contain(x, y)) continue;
		Point p = polygons[i].FindNear(x, y);
		double dx = p.x - x;
		double dy = p.y - y;
		double dxy = dx*dx + dy*dy;
		if (dis == -1 || dis > dxy)
		{
			nearpoint[0] = p.x;
			nearpoint[1] = p.y;
			dis = dxy;
		}
	}
	if (dis == -1) return 0;
	return nearpoint.data( );
}

//返回 0 说明无法两点之间无法直达，返回 1 说明可以直达
int32_t Path::CheckPath(double startx, double starty, double endx, double endy) {
	Point start = Point(startx, starty);
	Point end = Point(endx, endy);
	int32_t sIndex = -1;
	int32_t stIndex = -1;
	int32_t eIndex = -1;
	int32_t etIndex = -1;
	for (int32_t i = 0; i < polygons.size(); i++) {
		double minx = polygons[i].grid.minx;
		double maxx = polygons[i].grid.maxx;
		double miny = polygons[i].grid.miny;
		double maxy = polygons[i].grid.maxy;
		if ((startx >= minx &&startx <= maxx) && (starty >= miny &&starty <= maxy)) {
			sIndex = i;
			stIndex = polygons[sIndex].FindTriangle(start);
		}
		if ((endx >= minx && endx <= maxx) && (endy >= miny&&endy <= maxy)) {
			eIndex = i;
			etIndex = polygons[eIndex].FindTriangle(end);
		}
	}
	if (sIndex != eIndex || stIndex < 0 || etIndex < 0) return 0;


	double minx, maxx, miny, maxy;
	Polygon* p = &polygons[sIndex];
	minx = p->grid.minx;
	maxx = p->grid.maxx;
	miny = p->grid.miny;
	maxy = p->grid.maxy;
	int32_t gride = p->grid.gride;
	int32_t xnum = p->grid.xnum;
	int32_t ynum = p->grid.ynum;

	Hash edgehashs;

	if (start.x > end.x) {
		Point demo = start;
		start = end;
		end = demo;
	}
	int32_t xn1 = (int32_t)((start.x - minx) / gride);
	int32_t xn2 = (int32_t)((end.x - minx) / gride);
	int32_t yn1 = (int32_t)((start.y - miny) / gride);
	int32_t yn2 = (int32_t)((end.y - miny) / gride);
	if (xn1 == xn2) {
		if (yn1 > yn2) {
			yn1 = yn1 ^ yn2;
			yn2 = yn1 ^ yn2;
			yn1 = yn1 ^ yn2;
		}
		for (int32_t j = yn1; j <= yn2; j++)
		{
			if (j > ynum) break;
			int32_t edgepos = (xn1 >= xnum ? xnum - 1 : xn1) + (j >= ynum ? ynum - 1 : j)*xnum;
			if (edgepos < 0 || edgepos > xnum * ynum) continue;
			for (int32_t k = 0; k < p->grid.cells[edgepos].edges.size(); k++) {
				int32_t hashk = p->grid.cells[edgepos].edges[k];
				if (edgehashs.find(hashk) == edgehashs.end()) {
					edgehashs.insert(make_pair(hashk, edgepos));
				}
			}
		}
	}
	else {
		double y = start.y;
		double x = start.x;
		for (int32_t i = xn1; i <= xn2; i++)
		{
			double x3 = (i + 1) * gride + minx;
			if (x3 > end.x) x3 = end.x;
			double y3 = (start.x - x3)*(start.y - end.y) / (end.x - start.x) + start.y;

			int32_t cur_x = (int32_t)((x - minx) / gride);
			int32_t cur_y = (int32_t)((y - miny) / gride);
			int32_t next_y = (int32_t)((y3 - miny) / gride);
			if (cur_y > next_y) {
				cur_y = cur_y ^ next_y;
				next_y = cur_y ^ next_y;
				cur_y = cur_y ^ next_y;
			}
			for (int32_t j = cur_y; j <= next_y; j++)
			{
				if (j > ynum) break;
				int32_t edgepos = (cur_x >= xnum ? xnum - 1 : cur_x) + (j >= ynum ? ynum - 1 : j)*xnum;
				if (edgepos < 0 || edgepos > xnum * ynum) continue;
				for (int32_t k = 0; k < p->grid.cells[edgepos].edges.size(); k++) {
					int32_t hashk = p->grid.cells[edgepos].edges[k];
					if (edgehashs.find(hashk) == edgehashs.end()) {
						edgehashs.insert(make_pair(hashk, edgepos));
					}
				}
			}
			x = x3;
			y = y3;
		}
	}
	for (auto it = edgehashs.cbegin(); it != edgehashs.cend(); it++)
	{
		int32_t eposId = it->first;
		int32_t nextId = GetNextEposId(eposId, p);

		Point p1 = p->GetPoint(eposId);
		Point p2 = p->GetPoint(nextId);
		//TODO:求交点 start->face 与 p1->p2之间的交点
		if (Math::Meet(start, end, p1, p2)) {
			Point intep = Math::Inter(start, end, p1, p2);
			if ((intep.x == p1.x && intep.y == p1.y) || (intep.x == p2.x && intep.y == p2.y)
				|| (intep.x == start.x && intep.y == start.y) || (intep.x == end.x && intep.y == end.y)) {
				continue;
			}
			else {
				return 0;
			}
		}
	}
	
	return 1;
}

Path::~Path()
{
}
