

#include <iostream>
#include <Windows.h>


#include <vector>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;
typedef char*(*test)(char *);


typedef void*(*CreatePath)(double*, int);

typedef void(*SavePath)(void*, const char*);
typedef void*(*LoadPath)(const char*);

typedef void(*FreePath)(void*);
typedef const double*(*FindPath)(void*, double, double, double, double, int, int *);
typedef const double*(*FindCross)(void*, double, double, double, double);
typedef int(*CheckPath)(void*, double, double, double, double);
typedef const double*(*GetPoints)(void*, int*);
typedef const int*(*GetEdges)(void*, int*);


int main()
{
	HMODULE hModule = LoadLibrary("navmesh.dll");
	if (!hModule)
	{
		cout << "LoadLibrary ERROR!" << endl;
	}else
	{
		cout << "hModule Has Found!!!!" << endl;
	}
	

	ifstream infile;
	infile.open("./points.txt");
	if (!infile) {
		printf("unable to open the file!");
		exit(1);
	}
	printf("open the file!");
	string temp;
	int count = 0;
	double d;
	vector<double> cont;

	while (getline(infile, temp)) {
		if (temp.size() < 0) continue;
		d = atof(temp.c_str());
		cont.push_back(d);
		count++;
	}
	double* pos = new double[count];
	double* p = pos;

	for (int i = 0; i < cont.size(); i++) {
		*(double *)pos = cont[i];
		pos += 1;
	}

	CreatePath path = (CreatePath)GetProcAddress(hModule, "CreatePath");
	void *s1 = path(p, count);

	SavePath savep = (SavePath)GetProcAddress(hModule, "SavePath");
	savep(s1, "savepathdemo");

	LoadPath loadp = (LoadPath)GetProcAddress(hModule, "LoadPath");
	void * loadp1 = loadp("savepathdemo");
	

	FreeLibrary(hModule);

	system("Pause");
    return 0;
}

char* GetPointInfo(double pos[], char *p) {
	*(int *)p = sizeof(pos) / sizeof(double) / 2;
	p += sizeof(int);
	for (int i = 0; i < sizeof(pos) / sizeof(double); i++) {
		*(double *)p = pos[i];
		p += sizeof(double);
	}
	return p;
}