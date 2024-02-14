// NavPath.cpp : Defines the exported functions for the DLL application.
//

#include "stdio.h"
#include "Path.h"
#include "log.h"
using namespace std;

#if _WINDOWS
#define API extern "C" __declspec(dllexport)
#else
#define API extern "C"
#endif

API	char* echoString(char* cont) {
	return cont;
}


API  Path* CreatePath(double cont[], int32_t count) {
	int32_t totalcount = count;
	int32_t dvector = 0;		//cont的下标
	int32_t intnum = 0;			//所有int的个数
	int32_t needsubnum = 0;		//需要减去的int==0的个数
	//printf("create path...%d\n", totalcount);
	while (totalcount > 0) {
		int32_t childsize = (int32_t)cont[dvector++];
		
		int32_t parentlength = (int32_t)cont[dvector++];
		if (childsize != 0) {
			intnum += 2;
		}
		else {
			intnum += 1;
		}
		totalcount -= 2;
		dvector += parentlength * 2;
		totalcount -= parentlength * 2;
		if (childsize == 0) {
			intnum += 1;
		}
		else {
			needsubnum += childsize;
		}
		for (int32_t i = 0; i < childsize; i++) {
			int32_t childsizen = (int32_t)cont[dvector++];
			int32_t childlength = (int32_t)cont[dvector++];
			if (childsizen != 0) {
				intnum += 2;
			}
			else {
				intnum += 1;
			}
			totalcount -= 2;
			dvector += childlength * 2;
			totalcount -= childlength * 2;
		}
		//printf("childsizen %d %d\n", childsize, totalcount);
	}
	
	int32_t totalsize = count;		//double cont[] 的总长
	int32_t charlength = sizeof(double) * (count - intnum-needsubnum) + sizeof(int32_t)*intnum;
	char *pos = new char[charlength];
	char *p = pos;	//存储char的初始位置
	int32_t vector = 0;
	//printf("pos - p < charlength %d %d\n", pos - p, charlength);
	while (pos - p < charlength) {
		int32_t childn = (int32_t)cont[vector++];
		int32_t size = (int32_t)cont[vector++];

		*(int32_t *)pos = size;
		pos += sizeof(int32_t);
		
		for (int32_t i = vector; i < (size * 2 + vector); i++) {
			*(double *)pos = cont[i];
			pos += sizeof(double);
		}
		vector += size * 2;
		
		*(int32_t *)pos = childn;
		pos += sizeof(int32_t);

		if (childn > 0) {
			for (int32_t i = 0; i < childn; i++) {
				int32_t childsize = (int32_t)cont[++vector];
				*(int32_t *)pos = childsize;
				pos += sizeof(int32_t);
				for (int32_t j = vector+1; j < (childsize * 2 + vector+1); j++) {
					*(double *)pos = cont[j];
					pos += sizeof(double);
				}
				vector += childsize * 2;
				vector++;
			}
		}
		//printf("childn %d %d %d\n", childn, pos - p, charlength);
	}
	Path *path = new Path(p, charlength);
	delete[] p;

	return path;
} 

API  void SavePath(Path* handle, const char* resname) {
	handle->Save(resname);
}

API  void* LoadPath(const char* resname)
{
	LOGD("LoadPath %s", resname);
	Path* path = new Path();
	LOGD("LoadPath2 %s", resname);
	path->Load(resname);
	LOGD("LoadPath3 %s", resname);
	return path;
}

API  void FreePath(Path* handle) {
	delete handle;
}

API  const double* FindPath(Path* handle, double startx, double starty, double endx, double endy,int32_t far, int32_t* size) {
	Point start = Point(startx, starty);
	Point end = Point(endx, endy);
	const double* finalpath = handle->FindPaths(start, end, far==1, size);
	return finalpath;
}

API  const double* FindCross(Path* handle, double startx, double starty, double facex, double facey) {
	const double* crosspoint = handle->FindCross(startx, starty, facex, facey);
	return crosspoint;
}

API int32_t CheckPath(Path* handle, double startx, double starty, double endx, double endy) {
	int32_t result = handle->CheckPath(startx, starty, endx, endy);
	return result;
}

API const double* GetPoints(Path* path, int32_t* length) {
	const double* points = path->GetPoints(length);
	return points;
}

API const int32_t* GetSize(Path* path, int32_t* length) {
	const int32_t* points = path->GetSizes(length);
	return points;
}

API const int32_t* GetIndexs(Path* path, int32_t* length) {
	const int32_t* indexs = path->GetIndexs(length);
	return indexs;
}

API const double* FindNear(Path* handle, double x, double y)
{
	const double* point = handle->FindNear(x, y);
	return point;
}

API const int32_t DeleteFile( const char* filename ) {
	return remove(filename);
}

API const int32_t RenameFile(const char* oldname, const char* newname) {
	return rename(oldname, newname);
}