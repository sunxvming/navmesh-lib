#include "Polygon.h"
#include "log.h"
#include "Path.h"
//typedef std::tr1::unordered_map<int32_t, int32_t> Hash;
#define abs fabs
#define MAXDIS 99999.0f
#define USEDIJKSTRA true
Polygon::Polygon(double* pos, int32_t size)
{
	_ASSERT(size > 10);
	_ASSERT(size < MAXPOINT * 2 - 8);
	/*double minx = pos[0], miny = pos[1], maxx = pos[2], maxy = pos[3];
	points.push_back(Point(minx, miny));
	points.push_back(Point(maxx, miny));
	points.push_back(Point(maxx, maxy));
	points.push_back(Point(minx, maxy));*/
	for (int32_t ii = 0; ii < size; ii += 2)
	{
		points.push_back(Point(pos[ii], pos[ii + 1]));
	}
	Gen();
}

Polygon::Polygon()
{
}

Polygon Polygon::CreateFromShort(char* cont, int32_t*len)
{
	Polygon p;
	char* cont0 = cont;

	int32_t size = *(int32_t*)cont;
	cont += sizeof(int32_t);
	for (int32_t i = 0; i < size * 2; i++)
	{
		p.points.push_back(Point(((double *)cont)[i], ((double *)cont)[i + 1]));
		i++;
	}
	cont += sizeof(double)*size * 2;
	p.pointsnum.push_back(size);

	int32_t childn = *(int32_t*)cont;
	cont += sizeof(int32_t);
	//printf("create from short ==222222=====%d\n", childn);
	for (int32_t i = 0; i < childn; i++)
	{	
		int32_t childsize = *(int32_t*)cont;
		cont += sizeof(int32_t);
		//printf("create from short ==333333333=====%d\n", childsize);
		for (int32_t j = 0; j < childsize * 2; j++)
		{
			p.points.push_back(Point(((double *)cont)[j], ((double *)cont)[j + 1]));
			j++;
		}
		cont += sizeof(double) * childsize * 2;
		p.pointsnum.push_back(childsize);
	}
	*len = cont - cont0;
	p.Gen();
	return p;
}

void Polygon::Gen()
{
	//printf("start Delaunay\n");
	Delaunay();
	//printf("start GenExtData\n");
	GenExtData();
}
void Polygon::InitDijCache(FILE *file)
{
	int32_t size = points.size();
	vector<vector<DijCache>> dij1;
	vector<vector<int>> flag;
	dij1.resize(size);
	flag.resize(size);
	for (int32_t i = 0; i < size; i++)
	{
		dij1[i].resize(size);
		flag[i].resize(size);
		Point pi = points[i];
		for (int32_t j = 0; j < size; j++)
		{
			if (i == j)
			{
				dij1[i][j].preindex = i;
				dij1[i][j].dis = 0;
				continue;
			}
			flag[i][j] = 0;
			Point pj = points[j];
			dij1[i][j].dis = MAXDIS;
			//if ((i == 243 && j == 13) || (j == 243 && i == 13))
				//printf("check 240--- %f %f %f %f %d\n", pi.x, pi.y, pj.x, pj.y, FindTriangle(Point((pi.x + pj.x) / 2, (pi.y + pj.y) / 2)));
			int32_t ckpath = CheckPath(pi.x, pi.y, pj.x, pj.y);
			if (i != j && ((ckpath == 1 && FindTriangle(Point((pi.x + pj.x) / 2, (pi.y + pj.y) / 2)) >= 0) || ckpath == 2))
			{
				dij1[i][j].preindex = j;
				dij1[i][j].dis = (float)sqrt((pi.x - pj.x) * (pi.x - pj.x) + (pi.y - pj.y) * (pi.y - pj.y));
				//if ((i == 243 && j == 13) || (j == 243 && i == 13))
					//printf("dis %f %d %d\n", dij1[i][j].dis, i, j);
			}
		}
	}
	for (int32_t i = 0; i < size; i++)
	{
		float min;
		for (int32_t m = 1; m < size; m++)
		{
			min = MAXDIS;
			int32_t k = -1;
			for (int32_t j = 0; j < size; j++)//寻找当前最小的路径，中最小的权的定点
			{
				if (i == j)
					continue;
				if (flag[i][j] == 0 && dij1[i][j].dis < min)
				{
					min = dij1[i][j].dis;
					k = j;
				}
			}
			if (k == -1)
				break;
			flag[i][k] = 1;// 标记"顶点k"为已经获取到最短路径
			float temp;
			for (int32_t j = 0; j < size; j++)//修正当前最短路径和前驱顶点，也就是拿到顶点k的最短路径之后，更新未获取最短路径的顶点的最短路径和前驱顶点
			{
				if (i == j)
					continue;
				temp = min + dij1[k][j].dis;
				if (flag[i][j] == 0 && temp < dij1[i][j].dis)
				{
					//if ((i == 243 && j == 13) || (j == 243 && i == 13))
						//printf("temptemptemptemp %d %d %d %f \n", i, j, k, temp);
					dij1[i][j].dis = temp;
					dij1[i][j].preindex = dij1[i][k].preindex;
				}
			}
		}
	}
	for (int16_t i = 0; i < size; i++)
	{
		for (int16_t j = i + 1; j < size; j++)
		{
			//if ((i == 240 && j == 132) || (j == 240 && i == 132))
				//printf("check 240=== %d %d %d %f \n", i, j, dij1[i][j].preindex, dij1[i][j].dis);
			//if(i < 150)
			//printf("write init dij %d %d %d %f \n", i, j, dij1[i][j].preindex, dij1[i][j].dis);
			fwrite(&dij1[i][j].preindex, sizeof(int16_t), 1, file);
			//fwrite(&dij1[j][i].preindex, sizeof(int16_t), 1, file);
			fwrite(&dij1[i][j].dis, sizeof(float), 1, file);
		}
	}
	size = triangles.size();
	vector<vector<DijCache>> center1;
	vector<vector<int>> flagc;
	center1.resize(size);
	flagc.resize(size);
	for (int32_t i = 0; i < size; i++)
	{
		center1[i].resize(size);
		flagc[i].resize(size);
		//memset(flagc[i], 0, size * sizeof(int));
		Point pi = triangles[i].icenter;
		for (int32_t j = 0; j < size; j++)
		{
			if (i == j)
			{
				center1[i][j].preindex = i;
				center1[i][j].dis = 0;
				continue;
			}
			flagc[i][j] = 0;
			Point pj = triangles[j].icenter;
			center1[i][j].dis = MAXDIS;
			center1[i][j].preindex = -1;
			int32_t ckpath = CheckPath(pi.x, pi.y, pj.x, pj.y);
			if (i != j && ((ckpath == 1 && FindTriangle(Point((pi.x + pj.x) / 2, (pi.y + pj.y) / 2)) >= 0) || ckpath == 2))
			{
				center1[i][j].preindex = j;
				center1[i][j].dis = (float)sqrt((pi.x - pj.x) * (pi.x - pj.x) + (pi.y - pj.y) * (pi.y - pj.y));
			}
		}
	}
	/*for (int32_t i = 0; i < size; i++)
	{
		for (int32_t j = 0; j < size; j++)
		{
			printf("%d %d %d %f %d %d %f %d \n", i, j, center[i][j].preindex, center[i][j].dis, flagc[i][j]);
		}
		printf("\n");
	}*/
	for (int32_t i = 0; i < size; i++)
	{
		float min;
		for (int32_t m = 1; m < size; m++)
		{
			min = MAXDIS;
			//printf("before kkkkkkk %d %d \n", i, m);
			int16_t k = -1;
			for (int32_t j = 0; j < size; j++)//寻找当前最小的路径
			{
				if (i == j)
					continue;
				if (flagc[i][j] == 0 && center1[i][j].dis < min)
				{
					min = center1[i][j].dis;
					k = j;
				}
			}
			//printf("after kkkkkkk %f %d %d \n", min, i, k);
			if (k == -1)
				break;
			flagc[i][k] = 1;// 标记"顶点k"为已经获取到最短路径
			double temp;
			for (int32_t j = 0; j < size; j++)//修正当前最短路径和前驱顶点，也就是拿到顶点k的最短路径之后，更新未获取最短路径的顶点的最短路径和前驱顶点
			{
				if (i == j)
					continue;
				temp = min + center1[k][j].dis;
				if (flagc[i][j] == 0 && temp < center1[i][j].dis)
				{
					center1[i][j].dis = temp;
					center1[i][j].preindex = center1[i][k].preindex;
				}
			}
		}
	}
	for (int16_t i = 0; i < size; i++)
	{
		for (int16_t j = i + 1; j < size; j++)
		{
			//printf("write init dij %d %d %d %f \n", i, j, dij1[i][j].preindex, dij1[i][j].dis);
			fwrite(&center1[i][j].preindex, sizeof(int16_t), 1, file);
			//fwrite(&center1[j][i].preindex, sizeof(int16_t), 1, file);
			fwrite(&center1[i][j].dis, sizeof(float), 1, file);
		}
	}
}

void Polygon::CreateTriangle(Hash* eindexs, int32_t p1, int32_t p2, int32_t p3)
{
	triangles.push_back(Triangle(p1, p2, p3));
	int32_t triangle = triangles.size() - 1;
	triangles[triangle].edges[0] = CreateEdge(eindexs, triangle, p1, p2);
	triangles[triangle].edges[1] = CreateEdge(eindexs, triangle, p3, p2);
	triangles[triangle].edges[2] = CreateEdge(eindexs, triangle, p1, p3);
}

int32_t Polygon::CreateEdge(Hash* eindexs, int32_t triangle, int32_t p1, int32_t p2)
{
	int32_t k = PIndex(p1, p2);
	if (eindexs->find(k) == eindexs->end())
	{
		int32_t v = edges.size();
		edges.push_back(Edge(triangle, -1, p1, p2));
		eindexs->insert(make_pair(k, v));
		return v;
	}
	int32_t v = (*eindexs)[k];
	int32_t t1 = edges[v].triangles[0];
	if (t1 == v) return v;
	int32_t t2 = edges[v].triangles[1];
	if (t2 == v)return v;
	_ASSERT(t2 < 0);
	edges[v].triangles[1] = triangle;
	return v;
}


static inline int32_t visible(Point p1, Point p2, Point p3)			//p3是否在p1p2的右侧
{
	double x1 = p1.x, y1 = p1.y;
	double x2 = p2.x, y2 = p2.y;
	double x3 = p3.x, y3 = p3.y;
	return (x1 - x3)*(y2 - y3) - (y1 - y3)*(x2 - x3) > 0;
}

static inline int32_t visible(Point p2, Point p3, bool equal = true)		//p2是否在p1p3的右侧,equal，表示共线算不算/p2是否
{
	double x2 = p2.x, y2 = p2.y;
	double x3 = p3.x, y3 = p3.y;
	if (equal) return y3*x2 >= x3*y2;
	return y3*x2 > x3*y2;
}

static  inline double Angle(Point cen, Point first, Point second)
{
	double dx1, dx2, dy1, dy2;

	dx1 = first.x - cen.x;
	dy1 = first.y - cen.y;

	dx2 = second.x - cen.x;

	dy2 = second.y - cen.y;

	double c = (double)sqrt(dx1 * dx1 + dy1 * dy1) * (double)sqrt(dx2 * dx2 + dy2 * dy2);

	if (c == 0) return 0;

	return (dx1 * dx2 + dy1 * dy2) / c;
}

static inline int32_t lerp(int32_t v, int32_t min, int32_t max)
{
	if (v < min)
		v = min;
	if (v > max)
		v = max;
	return v;
}

int32_t Polygon::FindDT(Grid* grid, int32_t p1, int32_t p2)
{
	double x1 = points[p1].x, y1 = points[p1].y;
	double x2 = points[p2].x, y2 = points[p2].y;
	double x = (x1 + x2) / 2, y = (y1 + y2) / 2;

	double minx = grid->minx, miny = grid->miny;
	double gride = grid->gride;
	int32_t gx = (int32_t)((x - minx) / gride);
	int32_t gy = (int32_t)((y - miny) / gride);
	int32_t xnum = grid->xnum, ynum = grid->ynum;
	int32_t d = 0;
	int32_t p3 = -1;
	double angle3 = -1;
	Point point1 = GetPoint(p1), point2 = GetPoint(p2);
	while (1)
	{
		//printf("point1point1point1 %d \n", d);
		for (int32_t i = -d; i <= d; i++)

			for (int32_t j = -d; j <= d; j++)
			{
				int32_t pos = lerp(gx + i, 0, xnum - 1) + lerp(gy + j, 0, ynum - 1)*xnum;
				if (pos >= 0 && pos < (int32_t)grid->cells.size())
				{
					Cell c = grid->cells[pos];
					for (unsigned k = 0; k < c.points.size(); k++)
					{
						int32_t p = c.points[k];
						Point point = GetPoint(p);
						if (p1 != p && p2 != p &&visible(point1, point2, point))
						{
							bool flag = JudgeIsVisible(p1, p) && JudgeIsVisible(p2, p);
							if (flag) {
								double angle = Angle(point, point1, point2);
								if (p3 == -1 || angle < angle3)
								{
									angle3 = angle;
									p3 = p;
								}
							}

						}
					}
				}
			}
		//判断是否应该结束当前趟
		//printf("p3p3p3p3p3p3p3p3p3p3 %d \n", p3);
		if (p3 != -1)
		{
			Circle c(point1, point2, GetPoint(p3));
			Point cc = c.GetCenter();
			double radius = c.GetR();
			double l = cc.x - radius, r = cc.x + radius, t = cc.y - radius, b = cc.y + radius;
			int32_t lx = lerp((int32_t)((l - minx) / gride), 0, xnum - 1);
			int32_t rx = lerp((int32_t)((r - minx) / gride), 0, xnum - 1);
			int32_t ty = lerp((int32_t)((t - miny) / gride), 0, ynum - 1);
			int32_t by = lerp((int32_t)((b - miny) / gride), 0, ynum - 1);
			//printf("circleeeeeeeeeeeeeeeeee %d %d %d %d %d\n", lx, rx, ty, by, d);
			if ((gx - d) <= lx && (gx + d) >= rx && (gy - d) <= ty && (gy + d) >= by)
				break;
		}
		//_ASSERT(<);
		d++;
	}
	_ASSERT(p3 != -1);
	return p3;
}

Point Polygon::GetPoint(int32_t p)
{
	return points[p];
}

bool Polygon::IsIntersect(int32_t edgepos, int32_t pa1, int32_t p1) {
	Point pa = GetPoint(pa1);
	Point p = GetPoint(p1);
	Cell c = grid.cells[edgepos];
	bool flag = false;
	for (unsigned m = 0; m < c.edges.size(); m++) {
		int32_t eposId = c.edges[m];
		int32_t next_eposId = 0;
		int32_t sum = 0;
		for (int32_t i = 0; i < pointsnum.size(); i++) {
			sum += pointsnum[i];

			if (eposId == sum - 1) {
				next_eposId = eposId - (pointsnum[i] - 1);
				break;
			}
			else if (eposId < sum - 1) {
				next_eposId = eposId + 1;
				break;
			}
			else if (eposId > sum - 1) {
				continue;
			}
		}

		flag = Math::Meet(GetPoint(next_eposId), GetPoint(eposId), pa, p);
		if ((eposId == pa1) || (eposId == p1) || (next_eposId == pa1) || (next_eposId == p1))
		{
			flag = false;
		}
		if (flag)
			return false; 
	}
	return true;
}

//false代表pa 与p不可见，true代表可见
bool Polygon::JudgeIsVisible(int32_t pa1, int32_t p1) {
	Point pa = GetPoint(pa1);
	Point p = GetPoint(p1);
	//Hash passeline;		//存储两个点之间经过的所有的约束边
	//检查最大最小值合法性
	Point p0 = points[0];
	double minx = p0.x, miny = p0.y, maxx = p0.x, maxy = p0.y;
	for (auto it = points.cbegin() + 1; it != points.cend(); it++)
	{
		double x = it->x;
		double y = it->y;
		if (x < minx) minx = x;
		if (x > maxx)maxx = x;
		if (y < miny) miny = y;
		if (y > maxy)maxy = y;
	}

	double dx = maxx - minx, dy = maxy - miny;
	int32_t gride = (int32_t)sqrt(dx*dy / (points.size()));
	int32_t xnum = (int32_t)ceil(dx / gride);
	int32_t ynum = (int32_t)ceil(dy / gride);

	if (p.x > pa.x) {
		Point demo = p;
		p = pa;
		pa = demo;
	}
	int32_t xn1 = (int32_t)((p.x - minx) / gride);
	int32_t xn2 = (int32_t)((pa.x - minx) / gride);
	int32_t yn1 = (int32_t)((p.y - miny) / gride);
	int32_t yn2 = (int32_t)((pa.y - miny) / gride);

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
			if (!IsIntersect(edgepos, pa1, p1))
				return false;
		}
	}
	else {
		double y = p.y;
		double x = p.x;
		for (int32_t i = xn1; i <= xn2; i++)
		{
			double x3 = (i + 1) * gride + minx;
			if (x3 > pa.x) x3 = pa.x;
			double y3 = (p.x - x3)*(p.y - pa.y) / (pa.x - p.x) + p.y;

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
				if (!IsIntersect(edgepos, pa1, p1))
					return false;
			}
			x = x3;
			y = y3;
		}
	}
	return true;
}

static inline vector<int32_t>::iterator find(vector<int32_t>::iterator begin, vector<int32_t>::iterator end, int32_t v)
{
	while (begin != end)
	{
		if (*begin == v) return begin;
		begin++;
	}
	return end;
}

void Polygon::Delaunay()
{
	//检查最大最小值合法性
	Point p0 = points[0];
	double minx = p0.x, miny = p0.y, maxx = p0.x, maxy = p0.y;
	for (auto it = points.cbegin() + 1; it != points.cend(); it++)
	{
		double x = it->x;
		double y = it->y;
		if (x < minx) minx = x;
		if (x > maxx)maxx = x;
		if (y < miny) miny = y;
		if (y > maxy)maxy = y;
	}
	//printf("minxminxminxminx %f %f %f %f\n", minx, maxx, miny, maxy);
	double dx = maxx - minx, dy = maxy - miny;
	int32_t gride = (int32_t)sqrt(dx*dy / (points.size()));
	if (gride < 1)
		gride = 1;
	int32_t xnum = (int32_t)ceil(dx / gride);
	int32_t ynum = (int32_t)ceil(dy / gride);
	vector<Cell> cells(xnum*ynum);
	//printf("points.size()222222222222points.size() %f %d\n", dx*dy, points.size());
	for (auto it = points.cbegin(); it != points.cend(); it++)
	{

		double x = it->x;
		double y = it->y;

		int32_t xn = (int32_t)((x - minx) / gride);
		int32_t yn = (int32_t)(((y - miny) / gride));
		int32_t pos = (xn >= xnum ? xnum - 1 : xn) + (yn >= ynum ? ynum - 1 : yn)*xnum;
		//printf("push_back %d %d %d %d %f %f\n", gride, pos, xn, yn, x, y);
		int32_t point = it - points.cbegin();
		cells[pos].points.push_back(point);
		Point p = GetPoint(point);
		int32_t p1num = 0;
		int32_t sum = 0;
		for (int32_t i = 0; i < pointsnum.size(); i++) {
			
			sum += pointsnum[i];
			if(point == sum-1){
				p1num = point - (pointsnum[i] - 1);
				break;
			}
			else if(point < sum -1){
				p1num = point + 1;
				break;
			}
			else if (point > sum - 1) {
				continue;
			}
		}
		
		Point p1 = GetPoint(p1num);
		if (p.x > p1.x) {
			Point demo = p;
			p = p1;
			p1 = demo;
		}
		int32_t xn1 = (int32_t)((p.x - minx) / gride);
		int32_t xn2 = (int32_t)((p1.x - minx) / gride);
		int32_t yn1 = (int32_t)((p.y - miny) / gride);
		int32_t yn2 = (int32_t)((p1.y - miny) / gride);
		//if (xn1 == 12 || xn2 == 12 || yn1 == 6 || yn2 == 6)
	  /*if(p1.x > 161 && p1.x < 162 && p1.y > 92.5 && p1.y < 94)
			printf("pointpointpoint %d %d %d %d %f %f %f %f %d \n", xn1, xn2, yn1, yn2, p1.x, p1.y, p.x, p.y, point);
		if (p1.x > 169 && p1.x < 170 && p1.y > 94 && p1.y < 95)
			printf("p22ointpointpoint %d %d %d %d %f %f %f %f %d \n", xn1, xn2, yn1, yn2, p1.x, p1.y, p.x, p.y, point);
		if (p.x > 161 && p.x < 162 && p.y > 92.5 && p.y < 94)
			printf("p333ointpointpoint %d %d %d %d %f %f %f %f %d \n", xn1, xn2, yn1, yn2, p1.x, p1.y, p.x, p.y, point);
		if (p.x > 169 && p.x < 170 && p.y > 94 && p.y < 95)
			printf("p444ointpointpoint %d %d %d %d %f %f %f %f %d \n", xn1, xn2, yn1, yn2, p1.x, p1.y, p.x, p.y, point);*/
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
				cells[edgepos].edges.push_back(point);
			}
		}
		else {
			//将x值偏小的 y值保存起来
			double y = p.y;
			double x = p.x;
			int32_t cur_x = (int32_t)((x - minx) / gride);
			int32_t cur_y = (int32_t)((y - miny) / gride);
			for (int32_t i = xn1; i <= xn2; i++)
			{
				double x3 = (i + 1) * gride + minx;
				if (x3 > p1.x) x3 = p1.x;
				double y3 = (p.x - x3)*(p.y - p1.y) / (p1.x - p.x) + p.y;
				int32_t next_y = (int32_t)((y3 - miny) / gride);
				int32_t oldnexty = next_y;
				if (cur_y > next_y) {
					cur_y = cur_y ^ next_y;
					next_y = cur_y ^ next_y;
					cur_y = cur_y ^ next_y;
				}
				for (int32_t j = cur_y; j <= next_y; j++)
				{
					if (j > ynum) break;
					int32_t edgepos = (cur_x >= xnum ? xnum - 1 : cur_x) + (j >= ynum ? ynum - 1 : j)*xnum;
					cells[edgepos].edges.push_back(point);
				}
				x = x3;
				y = y3;
				cur_x++;
				cur_y = oldnexty;
			}
		}
		
	}

	grid.cells = cells;
	grid.gride = gride;
	grid.minx = minx;
	grid.maxx = maxx;
	grid.miny = miny;
	grid.maxy = maxy;
	grid.xnum = xnum;
	grid.ynum = ynum;

	Hash eindexs;
	Hash restrains;
	int32_t pos = 0;
	//printf("points.size()points.size() %d\n", points.size());
	for (unsigned i = 0; i < points.size(); i++)
	{
		int32_t num1 = i;
		int32_t sum = 0;
		int32_t num2 = 0;
		for (int32_t i = 0; i < pointsnum.size(); i++) {
			sum += pointsnum[i];

			if (num1 == sum - 1) {
				num2 = num1 - (pointsnum[i] - 1);
				break;
			}
			else if (num1 < sum - 1) {
				num2 = num1 + 1;
				break;
			}
			else if (num1 > sum - 1) {
				continue;
			}
		}
		restrains.insert(make_pair((int32_t)PIndex(num1, num2), 1));
	}
	vector<int32_t> es;
	int32_t p1 = 0, p2 = 1;
	int32_t e = PIndex(p1, p2);
	while (1)
	{
		//printf("while111111111\n");
		int32_t p3 = FindDT(&grid, p1, p2);
		//printf("while@@@@@@@ %d \n", p3);
		if (restrains.find((int32_t)PIndex(p1, p3)) == restrains.end())
		{
			vector<int32_t>::iterator it = find(es.begin(), es.end(), PIndex(p1, p3));
			if (it == es.end())
			{
				//printf("es.begin() %d\n", PIndex(p1, p3));
				es.push_back(PIndex(p1, p3));
			}
			else {
				//printf("es.erasees.erase==@@@@==\n");
				es.erase(it);
			}
		}
		//printf("while222222222222222\n");
		if (restrains.find((int32_t)PIndex(p2, p3)) == restrains.end())
		{
			vector<int32_t>::iterator it = find(es.begin(), es.end(), PIndex(p2, p3));
			if (it == es.end())
			{
				//printf("PIndex(p2, p3) %d\n", PIndex(p2, p3));
				es.push_back(PIndex(p2, p3));
			}
			else {
				es.erase(it);
			}
		}
		//printf("while3333333333333333333333\n");
		CreateTriangle(&eindexs, p1, p2, p3);
		//printf("es.empty() %d\n", es.size());
		if (es.empty()) { 
			break; 
		}
		e = *(es.end() - 1);
		es.pop_back();
		//SIndex(p1, p2, e);
		int32_t* points = edges[eindexs[e]].points;
		p1 = points[0], p2 = points[1];
	}
}

inline double Polygon::Distance(Point from, Point center, Point to)
{
	return distance(from, center) + distance(center, to);
}

double Polygon::Distance(Point from, int32_t edge, int32_t tindex, Point to)
{
	int32_t p0 = edges[edge].points[0];
	int32_t p1 = edges[edge].points[1];

	Point e0 = GetPoint(p0);
	Point e1 = GetPoint(p1);
	Point ecenter = Point((e0.x + e1.x) / 2, (e0.y + e1.y) / 2);

	int32_t p2 = triangles[tindex].p3;
	if (p2 == p0 || p2 == p1) p2 = triangles[tindex].p2;
	if (p2 == p0 || p2 == p1) p2 = triangles[tindex].p1;

	Point e2 = GetPoint(p2);
	Point ecenter2;
	if (distance(e0, to) < distance(e1, to))
	{
		ecenter2.x = (e0.x + e2.x) / 2;
		ecenter2.y = (e0.y + e2.y) / 2;
	}
	else {
		ecenter2.x = (e1.x + e2.x) / 2;
		ecenter2.y = (e1.y + e2.y) / 2;
	}
	return distance(from, ecenter) + distance(ecenter, ecenter2) + distance(ecenter2, to);
}

void Polygon::GenExtData()
{
	for (unsigned i = 0;  i < triangles.size(); i++)
	{
		Triangle t = triangles[i];
		t.GenExtData(this);
		triangles[i] = t;
	}

	for (unsigned i = 0; i < edges.size(); i++)
	{
		Edge e = edges[i];
		int32_t e0 = e.triangles[0];
		int32_t e1 = e.triangles[1];
		if (e0 >= 0 && e1 >= 0) {
			Triangle t0 = triangles[e0];
			Triangle t1 = triangles[e1];
			int32_t p0 = e.points[0], p1 = e.points[1];
			if (!visible(GetPoint(p0) - t0.icenter, GetPoint(p1) - t0.icenter, false))
			{
				e.points[0] = p1;
				e.points[1] = p0;
				edges[i] = e;
			}
		}
	}
	//printf("endnGenExtData\n");
}

vector<Line> Polygon::GetLines() {
	vector<Line> lines(edges.size());
	for (unsigned i = 0; i < edges.size(); i++)
	{
		Line line;
		line.p1 = GetPoint(edges[i].points[0]);
		line.p2 = GetPoint(edges[i].points[1]);
		if (edges[i].IsRestrain(this))
		{
			line.color[0] = 1.0;
			line.color[1] = 0.0;
			line.color[2] = 0.0;
		}
		else {
			line.color[0] = 0.0;
			line.color[1] = 1.0;
			line.color[2] = 0.0;
		}
		lines[i] = line;
	}
	return lines;
}

vector<Line> Polygon::GetGrideLines() {
	vector<Line> lines(grid.ynum + grid.xnum + 2);
	for (int32_t i = 0; i <= grid.xnum; i++)
	{
		Line line;
		Point pa1, pa2;
		double x = i*grid.gride + grid.minx;
		pa1.x = x;
		pa1.y = grid.miny;
		pa2.x = x;
		pa2.y = grid.miny + grid.ynum * grid.gride;
		line.p1 = pa1;
		line.p2 = pa2;

		line.color[0] = 0.2f;
		line.color[1] = 0.3f;
		line.color[2] = 0.3f;

		lines[i] = line;
	}

	for (int32_t j = 0; j <= grid.ynum; j++)
	{
		Line line;
		Point pa1, pa2;
		double y = j * grid.gride + grid.miny;
		pa1.x = grid.minx;
		pa1.y = y;
		pa2.x = grid.minx + grid.xnum * grid.gride;
		pa2.y = y;
		line.p1 = pa1;
		line.p2 = pa2;

		line.color[0] = 0.2f;
		line.color[1] = 0.3f;
		line.color[2] = 0.3f;

		lines[j + grid.xnum + 1] = line;
	}
	return lines;
}

vector<Point> Polygon::GetCenters() {
	vector<Point> centers(triangles.size());
	for (unsigned i = 0; i < triangles.size(); i++)
	{
		centers[i] = triangles[i].icenter;
	}
	return centers;
}

int32_t Polygon::Contain(double x, double y)
{
	double minx, maxx, miny, maxy;
	minx = grid.minx;
	maxx = grid.maxx;
	miny = grid.miny;
	maxy = grid.maxy;
	if (x<minx || x>maxx || y<miny || y>maxy)  return 0;
	return 1;
}

int32_t Polygon::FindTriangle(Point p)
{
	for (unsigned i = 0; i < triangles.size(); i++)
	{
		if (triangles[i].Contain(this, p)) return i;
	}
	return -1;
}

int32_t Polygon::CheckPath(double startx, double starty, double endx, double endy)
{
	Point start = Point(startx, starty);
	Point end = Point(endx, endy);
	double minx, maxx, miny, maxy;
	/*Polygon* p = &polygons[sIndex];*/
	minx = grid.minx;
	maxx = grid.maxx;
	miny = grid.miny;
	maxy = grid.maxy;
	int32_t gride = grid.gride;
	int32_t xnum = grid.xnum;
	int32_t ynum = grid.ynum;
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
	//printf("xnnnnnnnnn %f %f %f %f %d %d %d %d %d\n", startx, starty, endx, endy, xn1, xn2, yn1, yn2, gride);
	if (xn1 == xn2) {
		if (yn1 > yn2) {
			yn1 = yn1 ^ yn2;
			yn2 = yn1 ^ yn2;
			yn1 = yn1 ^ yn2;
			/*int32_t temp = yn1;
			yn1 = yn2;
			yn2 = temp;*/
		}
		for (int32_t j = yn1; j <= yn2; j++)
		{
			if (j > ynum) break;
			int32_t edgepos = (xn1 >= xnum ? xnum - 1 : xn1) + (j >= ynum ? ynum - 1 : j)*xnum;
			if (edgepos < 0 || edgepos > xnum * ynum) continue;
			//if (yn1 == 6 && xn1 == 12)
				//printf("edgesedgesedges %d %d %d %d \n", xnum, ynum, edgepos, grid.cells[edgepos].edges.size());
			for (int32_t k = 0; k < grid.cells[edgepos].edges.size(); k++) {
				int32_t hashk = grid.cells[edgepos].edges[k];
				if (edgehashs.find(hashk) == edgehashs.end()) {
					edgehashs.insert(make_pair(hashk, edgepos));
				}
			}
		}
	}
	else {
		double y = start.y;
		double x = start.x;
		//printf("xn1, xn1 %d %d\n", xn1, xn2);
		for (int32_t i = xn1; i <= xn2; i++)
		{
			double x3 = (i + 1) * gride + minx;
			if (x3 > end.x) x3 = end.x;
			double y3 = (start.x - x3)*(start.y - end.y) / (end.x - start.x) + start.y;

			//int32_t cur_x = (int32_t)((x - minx) / gride);
			int32_t cur_y = (int32_t)((y - miny) / gride);
			int32_t next_y = (int32_t)((y3 - miny) / gride);
			//printf("curxxxxxxxxxx %f %f %f %f %f %d \n", x, minx, y, miny, y3, gride);
			if (cur_y > next_y) {
				cur_y = cur_y ^ next_y;
				next_y = cur_y ^ next_y;
				cur_y = cur_y ^ next_y;
				/*int32_t temp = cur_y;
				cur_y = next_y;
				next_y = temp;*/
			}
			//printf("cur_y, next_y %d %d %d %d %f %f %f %f %f %f\n", i, cur_y, next_y, ynum, minx, miny, x, y, x3, y3);
			for (int32_t j = cur_y; j <= next_y; j++)
			{
				if (j > ynum) break;
				//int32_t edgepos = (cur_x >= xnum ? xnum - 1 : cur_x) + (j >= ynum ? ynum - 1 : j)*xnum;
				int32_t edgepos = (i >= xnum ? xnum - 1 : i) + (j >= ynum ? ynum - 1 : j)*xnum;
				//printf("cur_yyyyyyyyy %d %d %d %d %d\n", j, cur_x, xnum, edgepos, xnum * ynum);
				if (edgepos < 0 || edgepos > xnum * ynum) continue;
				//printf("grid.cells[edgepos].edges.size() %d\n", grid.cells[edgepos].edges.size());
				for (int32_t k = 0; k < grid.cells[edgepos].edges.size(); k++) {
					int32_t hashk = grid.cells[edgepos].edges[k];
					//printf("find?????? %d %d %d %d\n", edgepos, k, hashk, edgehashs.find(hashk) == edgehashs.end());
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
		int32_t nextId = GetNextEposId(eposId);
		Point p1 = GetPoint(eposId);
		Point p2 = GetPoint(nextId);
		//if(start.x >= 85.0f && start.x <= 85.4f && start.y > 125.5f && start.y < 125.7f)
			//printf("edgenow %d %d %f %f %f %f %d\n", eposId, nextId, p1.x, p1.y, p2.x, p2.y, Math::Meet(start, end, p1, p2));
		//TODO:求交点 start->face 与 p1->p2之间的交点
		//如果就是这条线段的端点的话,认为是相通的
		if (((abs(p1.x - startx) < 0.01f && abs(p1.y - starty) < 0.01f) && (abs(p2.x - endx) < 0.01f && abs(p2.y - endy) < 0.01f)) || 
			((abs(p2.x - startx) < 0.01f && abs(p2.y - starty) < 0.01f) && (abs(p1.x - endx) < 0.01f && abs(p1.y - endy) < 0.01f)))
		{
			//printf("the same line..............%d %d\n", eposId, nextId);
			return 2;
		}
		//printf("is meet?????? %f %f %f %f %d\n", p1.x, p1.y, p2.x, p2.y, Math::Meet(start, end, p1, p2));
		if (Math::Meet(start, end, p1, p2)) {
			Point intep = Math::Inter(start, end, p1, p2);
			//if (start.x >= 85.0f && start.x <= 85.4f && start.y > 125.5f && start.y < 125.7f)
			//printf("startx %f %f %f %f %f %f\n", intep.x, startx, intep.y, starty, endx, endy);
			if ((abs(intep.x - startx) < 0.01f && abs(intep.y - starty) < 0.01f) ||(abs(intep.x - endx) < 0.01f && abs(intep.y - endy) < 0.01f)) {  //交点如果是端点的话只有目标点才不算
				//if (start.x >= 85.0f && start.x <= 85.4f && start.y > 125.5f && start.y < 125.7f)
				//printf("continue@@@@@ %f %f %f %f\n", intep.x - startx, intep.y - starty, intep.x - endx, intep.y - endy);
				continue;
			}
			else 
			{
				//if (start.x >= 85.0f && start.x <= 85.4f && start.y > 125.5f && start.y < 125.7f)
					//printf("return 09\n");
				return 0;
			}
		}
	}
	//printf("checkpath succcc\n");
	return 1;
}

int32_t Polygon::GetNextEposId(int32_t eposId) {
	int32_t sum = 0;
	int32_t nextId = 0;
	for (int32_t i = 0; i < pointsnum.size(); i++)
	{
		sum += pointsnum[i];
		if (eposId == sum - 1) {
			nextId = eposId - (pointsnum[i] - 1);
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

int16_t Polygon::GetPreIndex(int16_t from, int16_t end, bool iscenter)
{
	if (iscenter)
	{
		if (from <= end)
		{
			return center[from][end].preindex;
		}
		else
		{
			//printf("center %d %d \n", center[end][from].preindex, end);
			if (center[end][from].preindex == from)
				return end;
			return GetPreIndex(from, center[end][from].preindex, iscenter);
		}
	}
	else
	{
		if (from <= end)
		{
			//printf("dddddddd %d %d %d %f \n", from, end, dij[from][end].preindex, dij[from][end].dis);
			return dij[from][end].preindex;
		}
		else
		{
			//printf("dij %d %d %d \n", dij[end][from].preindex, from, end);
			if (dij[end][from].preindex == from)
				return end;
			return GetPreIndex(from, dij[end][from].preindex, iscenter);
		}
	}
}

float Polygon::GetDijDis(Point from, Point to, int16_t pointfrom, int16_t pointto)
{
	Point pto = GetPoint(pointto);
	if (CheckPath(from.x, from.y, pto.x, pto.y) >= 1)
	{
		pointfrom = pointto;
	}
	int16_t pretoindex = GetPreIndex(pointto, pointfrom, false);
	while (pretoindex != pointfrom)
	{
		Point pretopoint = GetPoint(pretoindex);
		if (CheckPath(from.x, from.y, pretopoint.x, pretopoint.y) >= 1)
		{
			pointfrom = pretoindex;
			break;
		}
		pretoindex = GetPreIndex(pretoindex, pointfrom, false);
	}
	Point pfrom = GetPoint(pointfrom);
	if (CheckPath(pfrom.x, pfrom.y, to.x, to.y) >= 1)
	{
		pointto = pointfrom;
	}
	pretoindex = GetPreIndex(pointfrom, pointto, false);
	while (pretoindex != pointto)
	{
		Point pretopoint = GetPoint(pretoindex);
		if (CheckPath(pretopoint.x, pretopoint.y, to.x, to.y) >= 1)
		{
			pointto = pretoindex;
			break;
		}
		pretoindex = GetPreIndex(pretoindex, pointto, false);
	}
	//ways.push_back(GetPoint(pointfrom));
	Point pointf = GetPoint(pointfrom);
	float disnow = (float)sqrt((from.x - pointf.x) * (from.x - pointf.x) + (from.y - pointf.y) * (from.y - pointf.y));
	//printf("push_back %d %f %f\n", pointfrom, GetPoint(pointfrom).x, GetPoint(pointfrom).y);
	disnow += (pointto >= pointfrom ? dij[pointfrom][pointto].dis : dij[pointto][pointfrom].dis);
	Point pointt = GetPoint(pointto);
	disnow += (float)sqrt((to.x - pointt.x) * (to.x - pointt.x) + (to.y - pointt.y) * (to.y - pointt.y));
	return disnow;
}

//true时为中点寻路，false为拐点寻路)
vector<Point> Polygon::FindPath(Point from, Point to, bool isturn)
{
	int32_t tfrom = FindTriangle(from);
	int32_t tto = FindTriangle(to);
	vector<Point> ways;
	if (tfrom < 0) return ways;
	if (tto < 0) return ways;
	ways.push_back(from);
	//printf("checkpath@@@@@@@@@ %d", CheckPath(from.x, from.y, to.x, to.y));
	if (tfrom == tto || CheckPath(from.x, from.y, to.x, to.y) >= 1)
	{
		ways.push_back(to);
		return ways;
	}
	if (USEDIJKSTRA) //DIJKSTRA
	{
		if (isturn) //中点寻路
		{
			//Point pointto = triangles[tto].icenter;
			//printf("start isturnnnn %d %d \n", tfrom, tto);
			int16_t pretoindex = GetPreIndex(tfrom, tto, true);
			while (true)
			{
				Point pretopoint = triangles[pretoindex].icenter;
				//printf("checkpath11111111111=== %d %d %d\n", tfrom, tto, pretoindex);
				if (CheckPath(from.x, from.y, pretopoint.x, pretopoint.y) >= 1)
				{
					tfrom = pretoindex;
					pretoindex = GetPreIndex(pretoindex, tto, true);
					if (pretoindex == tto)
					{
						pretopoint = triangles[pretoindex].icenter;
						if (CheckPath(from.x, from.y, pretopoint.x, pretopoint.y) >= 1)
							tfrom = tto;
						break;
					}
				}
				else
					break;
			}
			//Point pointfrom = triangles[tfrom].icenter;
			/*if (CheckPath(pointfrom.x, pointfrom.y, to.x, to.y) >= 1)
			{
				tto = tfrom;
			}*/
			pretoindex = GetPreIndex(tto, tfrom, true);
			while (true)
			{
				Point pretopoint = triangles[pretoindex].icenter;
				//printf("checkpath222222222222222222222 %d %d %d %d\n", pretoindex, tto, tfrom, CheckPath(pretopoint.x, pretopoint.y, to.x, to.y));
				if (CheckPath(pretopoint.x, pretopoint.y, to.x, to.y) >= 1)
				{
					tto = pretoindex;
					pretoindex = GetPreIndex(pretoindex, tfrom, true);
					if (pretoindex == tfrom)
					{
						pretopoint = triangles[pretoindex].icenter;
						//printf("cherckpaaaaaaa %d \n", pretoindex);
						if (CheckPath(to.x, to.y, pretopoint.x, pretopoint.y) >= 1)
							tto = pretoindex;
						break;
					}
				}
				else
					break;
			}
			//Point ptfrom = Point(triangles[tfrom].icenter.x, triangles[tfrom].icenter.y);
			//Point ptfromnext = triangles[GetPreIndex(tfrom, tto, true)].icenter;
			Point ptfrom = triangles[tfrom].icenter;
			if (tfrom != tto)
			{
				Point ptnext = triangles[GetPreIndex(tfrom, tto, true)].icenter;
				if (Math::pointmul(from, ptfrom, ptnext) < 0) //夹角大于90
				{
					Point ppoint = Math::pendalpoint(from, ptfrom, ptnext);
					double fpdis = max(4.0, 0.04 * ((ppoint.x - from.x) * (ppoint.x - from.x) + (ppoint.y - from.y) * (ppoint.y - from.y)));
					//printf("checkpath3333333333333333333333\n");
					if (CheckPath(from.x, from.y, ppoint.x, ppoint.y) >= 1 && CheckPath(ptnext.x, ptnext.y, ppoint.x, ppoint.y) >= 1)
					{
						Point result = Point((ppoint.x + ptnext.x) / 2, (ppoint.y + ptnext.y) / 2);
						bool judge = false;
						while ((result.x - ppoint.x)*(result.x - ppoint.x) + (result.y - ppoint.y)*(result.y - ppoint.y) >= fpdis)
						{
							//printf("1111111111111111111");
							//printf("checkpath44444444444444444444444\n");
							if (CheckPath(from.x, from.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ppoint.x / 2 + result.x / 2;
							result.y = ppoint.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptfrom = result;
						}
						else
						{
							ptfrom = ppoint;
						}
					}
					else
					{
						Point result = Point((ppoint.x + ptfrom.x) / 2, (ppoint.y + ptfrom.y) / 2);
						bool judge = false;
						while ((result.x - ptfrom.x)*(result.x - ptfrom.x) + (result.y - ptfrom.y)*(result.y - ptfrom.y) >= fpdis)
						{
							//printf("22222222222222222");
							//printf("checkpath55555555555555555555\n");
							if (CheckPath(from.x, from.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ptfrom.x / 2 + result.x / 2;
							result.y = ptfrom.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptfrom = result;
						}
					}
				}
				ways.push_back(ptfrom);
				int32_t tpre = tfrom;
				while (tpre != tto)
				{
					tpre = GetPreIndex(tpre, tto, true);
					if (tpre == tto)
						break;
					ways.push_back(triangles[tpre].icenter);
				}
				int32_t ttopre = GetPreIndex(tto, tfrom, true);
				if (ttopre == tfrom)
				{
					ptnext = ptfrom;
				}
				else
				{
					ptnext = triangles[ttopre].icenter;
				}
				Point ptend = triangles[tto].icenter;
				if (Math::pointmul(to, ptend, ptnext) < 0) //夹角大于90
				{
					Point ppoint = Math::pendalpoint(to, ptend, ptnext);
					double fpdis = max(4.0, 0.04 * ((ppoint.x - to.x) * (ppoint.x - to.x) + (ppoint.y - to.y) * (ppoint.y - to.y)));
					//printf("checkpath666666666666666666666666\n");
					if (CheckPath(to.x, to.y, ppoint.x, ppoint.y) >= 1 && CheckPath(ptnext.x, ptnext.y, ppoint.x, ppoint.y) >= 1)
					{
						Point result = Point((ppoint.x + ptnext.x) / 2, (ppoint.y + ptnext.y) / 2);
						bool judge = false;
						while ((result.x - ppoint.x)*(result.x - ppoint.x) + (result.y - ppoint.y)*(result.y - ppoint.y) >= fpdis)
						{
							//printf("3333333333333333333");
							//printf("checkpath7777777777777777777\n");
							if (CheckPath(to.x, to.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ppoint.x / 2 + result.x / 2;
							result.y = ppoint.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptend = result;
						}
						else
						{
							ptend = ppoint;
						}
					}
					else
					{
						Point result = Point((ppoint.x + ptend.x) / 2, (ppoint.y + ptend.y) / 2);
						bool judge = false;
						while ((result.x - ptend.x)*(result.x - ptend.x) + (result.y - ptend.y)*(result.y - ptend.y) >= fpdis)
						{
							//printf("4444444444444444444444");
							//printf("checkpath8888888888888888888\n");
							if (CheckPath(to.x, to.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ptend.x / 2 + result.x / 2;
							result.y = ptend.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptend = result;
						}
					}
				}
				ways.push_back(ptend);
				ways.push_back(to);
				return ways;
			}
			else
			{
				//printf("check from == to %f\n", Math::pointmul(from, ptfrom, to));
				if (Math::pointmul(from, ptfrom, to) < 0) //夹角大于90
				{
					Point ppoint = Math::pendalpoint(from, ptfrom, to);
					double fpdis = max(4.0, 0.04 * ((ppoint.x - from.x) * (ppoint.x - from.x) + (ppoint.y - from.y) * (ppoint.y - from.y)));
					//printf("checkpath9999999999999999999999 %f %f \n", ppoint.x, ppoint.y);
					if (CheckPath(from.x, from.y, ppoint.x, ppoint.y) >= 1 && CheckPath(to.x, to.y, ppoint.x, ppoint.y) >= 1)
					{
						Point result = Point((ppoint.x + to.x) / 2, (ppoint.y + to.y) / 2);
						bool judge = false;
						while ((result.x - ppoint.x)*(result.x - ppoint.x) + (result.y - ppoint.y)*(result.y - ppoint.y) >= fpdis)
						{
							//printf("checkpath111=== %d \n", CheckPath(from.x, from.y, result.x, result.y));
							//printf("55555555555555555");
							if (CheckPath(from.x, from.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ppoint.x / 2 + result.x / 2;
							result.y = ppoint.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptfrom = result;
						}
						else
						{
							ptfrom = ppoint;
						}
					}
					else
					{
						Point result = Point((ppoint.x + ptfrom.x) / 2, (ppoint.y + ptfrom.y) / 2);
						bool judge = false;
						while ((result.x - ptfrom.x)*(result.x - ptfrom.x) + (result.y - ptfrom.y)*(result.y - ptfrom.y) >= fpdis)
						{
							//printf("6666666666666666666");
							//printf("checkpath10101101010101010101010\n");
							if (CheckPath(from.x, from.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ptfrom.x / 2 + result.x / 2;
							result.y = ptfrom.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptfrom = result;
						}
					}
				}
				if (Math::pointmul(to, ptfrom, from) < 0) //夹角大于90
				{
					Point ppoint = Math::pendalpoint(to, ptfrom, from);
					double fpdis = max(4.0, 0.04 * ((ppoint.x - to.x) * (ppoint.x - to.x) + (ppoint.y - to.y) * (ppoint.y - to.y)));
					//printf("checkpath====2222222222222222222222\n");
					if (CheckPath(to.x, to.y, ppoint.x, ppoint.y) >= 1 && CheckPath(from.x, from.y, ppoint.x, ppoint.y) >= 1)
					{
						Point result = Point((ppoint.x + from.x) / 2, (ppoint.y + from.y) / 2);
						bool judge = false;
						while ((result.x - ppoint.x)*(result.x - ppoint.x) + (result.y - ppoint.y)*(result.y - ppoint.y) >= fpdis)
						{
							//printf("checkpath111=== %d \n", CheckPath(from.x, from.y, result.x, result.y));
							//printf("777777777777777777777");
							//printf("checkpath====3333333333333333333333\n");
							if (CheckPath(to.x, to.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ppoint.x / 2 + result.x / 2;
							result.y = ppoint.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptfrom = result;
						}
						else
						{
							ptfrom = ppoint;
						}
					}
					else
					{
						Point result = Point((ppoint.x + ptfrom.x) / 2, (ppoint.y + ptfrom.y) / 2);
						bool judge = false;
						while ((result.x - ptfrom.x)*(result.x - ptfrom.x) + (result.y - ptfrom.y)*(result.y - ptfrom.y) >= fpdis)
						{
							//printf("888888888888888888888888");
							//printf("checkpath====4444444444444444444444\n");
							if (CheckPath(to.x, to.y, result.x, result.y) >= 1)
							{
								judge = true;
								break;
							}
							result.x = ptfrom.x / 2 + result.x / 2;
							result.y = ptfrom.y / 2 + result.y / 2;
						}
						if (judge)
						{
							ptfrom = result;
						}
					}
				}
				ways.push_back(ptfrom);
				ways.push_back(to);
				return ways;
			}
		}
		else
		{
			Triangle trifrom = triangles[tfrom];
			Triangle trito = triangles[tto];
			int16_t pointfrom = -1, pointto = -1;
			float min = MAXDIS;
			for (int32_t i = 1; i <= 3; i++)
			{
				int16_t pf;
				if (i == 1)
					pf = trifrom.p1;
				else
					if (i == 2)
						pf = trifrom.p2;
					else
						pf = trifrom.p3;
				for (int32_t j = 1; j <= 3; j++)
				{
					int16_t pt;
					if (j == 1)
						pt = trito.p1;
					else
						if (j == 2)
							pt = trito.p2;
						else
							pt = trito.p3;
					Point pointf = points[pf];
					Point pointt = points[pt];
					/*float disnow = (float)sqrt((from.x - pointf.x) * (from.x - pointf.x) + (from.y - pointf.y) * (from.y - pointf.y)) +
									(float)sqrt((to.x - pointt.x) * (to.x - pointt.x) + (to.y - pointt.y) * (to.y - pointt.y));*/
					float disnow = GetDijDis(from, to, pf, pt);
					//disnow += (pt >= pf ? dij[pf][pt].dis : dij[pt][pf].dis);
					//printf("disnow...... %d %d %f %f\n", pf, pt, disnow, min);
					if (min > disnow)
					{
						min = disnow;
						pointfrom = pf;
						pointto = pt;
						//printf("dispolynowwwwwwwwww %f %d %d \n", min, pointfrom, pointto);
					}
				}
			}
			//printf("pointfrom %d \n", pointfrom);
			if (pointfrom == -1)
				return ways;
			//printf("pointfromm===after %d %d %d\n", pointfrom, pointto);
			Point pto = GetPoint(pointto);
			if (CheckPath(from.x, from.y, pto.x, pto.y) >= 1)
			{
				//printf("pointfromm===after %d %d \n", pointfrom, pointto);
				pointfrom = pointto;
			}
			int16_t pretoindex = GetPreIndex(pointto, pointfrom, false);
			while (pretoindex != pointfrom)
			{
				//printf("pointfrommmmmmmmmmm%d %d \n", pretoindex, pointfrom);
				Point pretopoint = GetPoint(pretoindex);
				if (CheckPath(from.x, from.y, pretopoint.x, pretopoint.y) >= 1)
				{
					//printf("pointfrommm333333333mmmmmmmm%d %d \n", pretoindex, pointfrom);
					pointfrom = pretoindex;
					break;
				}
				pretoindex = GetPreIndex(pretoindex, pointfrom, false);
			}
			Point pfrom = GetPoint(pointfrom);
			//printf("first pointto %d %f %f \n", pointfrom, pfrom.x, pfrom.y);
			if (CheckPath(pfrom.x, pfrom.y, to.x, to.y) >= 1)
			{
				//printf("pointfrommmm2222222222222mmmmmmm%d %d \n", pointto, pointfrom);
				pointto = pointfrom;
			}
			pretoindex = GetPreIndex(pointfrom, pointto, false);
			while (pretoindex != pointto)
			{
				Point pretopoint = GetPoint(pretoindex);
				//printf("checktooooooooo %d %d\n", pretoindex, pointto);
				if (CheckPath(pretopoint.x, pretopoint.y, to.x, to.y) >= 1)
				{
					//printf("checkto22222222222222oooooooo %d %d\n", pretoindex, pointto);
					pointto = pretoindex;
					break;
				}
				pretoindex = GetPreIndex(pretoindex, pointto, false);
			}
			ways.push_back(GetPoint(pointfrom));
			//printf("push_back %d %f %f\n", pointfrom, GetPoint(pointfrom).x, GetPoint(pointfrom).y);
			while (pointfrom != pointto)
			{
				pointfrom = GetPreIndex(pointfrom, pointto, false);
				ways.push_back(GetPoint(pointfrom));
				//printf("push_back@ %d %d %f %f\n", pointfrom, pointto, GetPoint(pointfrom).x, GetPoint(pointfrom).y);
			}
			//printf("push_back## %f %f\n", to.x, to.y);
			ways.push_back(to);
			return ways;
		}
	}
	else //A*
	{
		vector<int32_t> ts;
		ts.push_back(tfrom);
		vector<int32_t> es;
		int32_t start = tfrom;
		Hash visited;
		visited.insert(make_pair(tfrom, 1));
		Point startp = from;
		while (!ts.empty())
		{
			start = *(ts.cend() - 1);
			double weight = -1;
			int32_t next = -1;
			int32_t e;
			for (int32_t i = 0; i < 3; i++)
			{
				int32_t eindex = triangles[start].edges[i];
				Edge edge = edges[eindex];
				int32_t nextt = edge.triangles[0] == start ? edge.triangles[1] : edge.triangles[0];
				if (nextt >= 0)
				{
					if (nextt == tto)
					{
						next = tto;
						e = eindex;
						break;
					}
					if (visited.find(nextt) == visited.end())
					{
						Point e1 = GetPoint(edge.points[0]);
						Point e2 = GetPoint(edge.points[1]);
						Point ecenter = Point((e1.x + e2.x) / 2, (e1.y + e2.y) / 2);
						Point p = triangles[nextt].icenter;
						//double d = Distance(startp, eindex, nextt, to);
						double d = Distance(p, ecenter, to);
						if (next < 0) next = nextt, weight = d, startp = ecenter, e = eindex;
						if (weight > d) next = nextt, weight = d, startp = ecenter, e = eindex;
					}
				}
			}
			if (next == tto)
			{
				ts.push_back(next);
				es.push_back(e);
				break;
			}
			else if (next >= 0)
			{
				visited.insert(make_pair(next, 1));
				ts.push_back(next);
				es.push_back(e);
			}
			else {
				ts.pop_back();
				es.pop_back();
			}
		}

		if (ts.empty()) return ways;
		int32_t t = 0;
		int32_t size = es.size();
		if (isturn) {
			//中点寻路算法
			for (int32_t i = 0; i < es.size(); i++) {
				Edge e = edges[es[i]];
				Point p1 = GetPoint(e.points[0]);
				Point p2 = GetPoint(e.points[1]);
				Point cp = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
				ways.push_back(cp);
			}
			ways.push_back(to);
		}
		else {
			while (t < size)
			{
				Point p = *(ways.cend() - 1);
				Edge e = edges[es[t]];
				int32_t ep0, ep1;
				if (ts[t] == e.triangles[0])
				{
					ep0 = e.points[0], ep1 = e.points[1];
				}
				else {
					ep0 = e.points[1], ep1 = e.points[0];
				}
				Point p1 = GetPoint(ep0), p2 = GetPoint(ep1);
				Point line1 = p1 - p, line2 = p2 - p;
				int32_t i = t + 1;
				int32_t t1 = i, t2 = i;
				for (; i < size; i++)
				{
					Edge nxte = edges[es[i]];
					int32_t nep0, nep1;
					if (ts[i] == nxte.triangles[0])
					{
						nep0 = nxte.points[0], nep1 = nxte.points[1];
					}
					else {
						nep0 = nxte.points[1], nep1 = nxte.points[0];
					}
					Point np1 = GetPoint(nep0), np2 = GetPoint(nep1);
					if (!visible(line1, np2 - p, true)) {//np2是否在上拐点的左侧(np2是否在pep0左侧）
						ways.push_back(GetPoint(ep0));
						t1++;
						t = nep0 == ep0 ? (i++, t1 + 1) : t1;
						break;
					}
					else if (visible(line1, np1 - p, true))
					{
						line1 = np1 - p;
						ep0 = nep0;
						t1 = i;
					}

					if (visible(line2, np1 - p, false)) { //np1是否在下拐点的右侧
						ways.push_back(GetPoint(ep1));
						t2++;
						t = nep1 == ep1 ? (i++, t2 + 1) : t2;
						break;
					}
					else if (!visible(line2, np2 - p, false))
					{
						line2 = np2 - p;
						ep1 = nep1;
						t2 = i;
					}
				}

				if (i >= size - 1)
				{
					if (!visible(line1, to - p, true))
					{
						ways.push_back(GetPoint(ep0));
						t = t1;
					}
					else if (visible(line2, to - p, false))
					{
						ways.push_back(GetPoint(ep1));
						t = t2;
					}
					if (t >= size - 1)
					{
						ways.push_back(to);

					}
					t++;
				}

			}
		}//End else
		return ways;
	}
}

Point Polygon::FindNear(double x, double y)
{
	//Point vector = Point(facex - startx, facey - starty);		//向量
	double minx, maxx, miny, maxy;
	minx = grid.minx;
	maxx = grid.maxx;
	miny = grid.miny;
	maxy = grid.maxy;
	int32_t gride = grid.gride;
	int32_t xnum = grid.xnum;
	int32_t ynum = grid.ynum;
	Point pt(x, y);
	int32_t gx = (int32_t)((x - minx) / gride);
	int32_t gy = (int32_t)((y- miny) / gride);
	int32_t  d = 0;
	int32_t maxd = max(xnum, ynum);
	tmpedges.clear();
	tmpcells.clear();
	int32_t edge = -1;
	double mindis = -1;
	while (1)
	{
		for (int32_t i = -d; i <= d; i++)
			for (int32_t j = -d; j <= d; j++)
			{
				int32_t pos = lerp(gx + i, 0, xnum - 1) + lerp(gy + j, 0, ynum - 1)*xnum;
				if (tmpcells.find(pos) != tmpcells.end())
					break;
				if (pos >= 0 && pos < (int32_t)grid.cells.size())
				{
					Cell c = grid.cells[pos];
					for (unsigned k = 0; k < c.edges.size(); k++)
					{
						int32_t e = c.edges[k];
						if (tmpedges.find(e) != tmpedges.end())
							continue;
						double dis = edges[e].Distance(this, pt);
						if (edge == -1 || mindis > dis)
						{
							mindis = dis;
							edge = e;
							maxd = (int32_t)(mindis / gride) + 2;
						}
						tmpedges.insert(make_pair(e, 1));
					}
				}
				tmpcells.insert(make_pair(pos, 1));
			}
		if (gx - d <= 0 && gx + d >= xnum && gy - d <= 0 && gy + d >= ynum)
			break;
		if (d >= maxd)
			break;
		d++;
	}
	_ASSERT(edge != -1);
	return edges[edge].FindPoint(this, pt);
}

void Polygon::Save(FILE* file) {
	int32_t pointsLen = points.size();
	fwrite(&(pointsLen), sizeof(int32_t), 1, file);
	fwrite(points.data(), sizeof(Point), pointsLen, file);
	//printf("pointslen %d\n", pointsLen);
	int32_t trianglesLen = triangles.size();
	fwrite(&(trianglesLen), sizeof(int32_t), 1, file);
	fwrite(triangles.data(), sizeof(Triangle), trianglesLen, file);
	//printf("trianglesLen %d\n", trianglesLen);
	int32_t edgeLen = edges.size();
	fwrite(&(edgeLen), sizeof(int32_t), 1, file);
	fwrite(edges.data(), sizeof(Edge), edgeLen, file);
	//printf("edgeLen %d\n", edgeLen);
	int32_t pLen = pointsnum.size();
	fwrite(&(pLen), sizeof(int32_t), 1, file);
	fwrite(pointsnum.data(), sizeof(int32_t), pLen, file);
	//printf("pLen %d\n", pLen);
	int32_t cellsize = grid.xnum * grid.ynum;
	//printf("cellsize %d\n", cellsize);
	fwrite(&(cellsize), sizeof(int32_t), 1, file);
	for (int32_t i = 0; i < cellsize; i++) {
		int32_t poLen = grid.cells[i].points.size();
		fwrite(&(poLen), sizeof(int32_t), 1, file);
		fwrite(grid.cells[i].points.data(), sizeof(int32_t), poLen, file);

 		int32_t edLen = grid.cells[i].edges.size();
		fwrite(&(edLen), sizeof(int32_t), 1, file);
		fwrite(grid.cells[i].edges.data(), sizeof(int32_t), edLen, file);
	}
	fwrite(&(grid.gride), sizeof(int32_t), 1, file);
	fwrite(&(grid.minx), sizeof(double), 1, file);
	fwrite(&(grid.maxx), sizeof(double), 1, file);
	fwrite(&(grid.miny), sizeof(double), 1, file);
	fwrite(&(grid.maxy), sizeof(double), 1, file);
	fwrite(&(grid.xnum), sizeof(int32_t), 1, file);
	fwrite(&(grid.ynum), sizeof(int32_t), 1, file);
	//printf("write grid=======%f %f %f %f \n", grid.minx, grid.maxx, grid.miny, grid.maxy);
	InitDijCache(file);
}

Polygon::~Polygon()
{
}
