#ifndef WIRE_H
#define WIRE_H
#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include "GL/glut.h"
#include "GL/glu.h"
#include "GL/glut.h"
//#include "Wire.h"
#include <math.h>
#include <time.h>
//#include "susan.h"

#define Examplifier 10000

using namespace std;



struct Vertex{
	double x,y,z;
	double nx,ny,nz;
};

struct Face{
	int v1,v2,v3;
	double nx,ny,nz;
};


struct Point{
	double x,y;
};

struct Point3{
	double x,y,z;
};

struct Point4{
	int i;
	int x,y;
	//cost z;
	int z;
/*	Point4* operator=(const Point4 &rhs)
	{
		this->i=rhs.i;
		this->x=rhs.x;
		this->y=rhs.y;
		this->z=rhs.z;
		return this;
	}*/

};

struct Vector{
	double x,y,w;
	};

typedef float Matrix[16];

/*
Vector mult(Matrix m, Vector &v)
{
	Vector result;
    result.x =m[0]*v.x+m[4]*v.y+m[8]*v.z+m[12]*v.w;
    result.y =m[1]*v.x+m[5]*v.y+m[9]*v.z+m[13]*v.w;
    result.z =m[2]*v.x+m[6]*v.y+m[10]*v.z+m[14]*v.w;
    result.w =m[3]*v.x+m[7]*v.y+m[11]*v.z+m[15]*v.w;
	return result;
}

inline
void dumpMatrix(Matrix m)
{
    for(int i = 0; i < 4; ++i)
    {
        fprintf(stderr, "%f %f %f %f\n",
                m[i * 4 + 0],
                m[i * 4 + 1],
                m[i * 4 + 2],
                m[i * 4 + 3]);
    }
}
*/
//Constant Macros
#define PI 3.1415926

//extern vector<Point> cubic_spline_interp(vector<Point> q,vector<double> u);
//extern vector<double> chordparametrize( vector <Point> q);
//extern vector<Point>  uniformsample(vector<Point> ctrlp,vector<double> uP);
extern void disp(vector<Point> T,int window);
extern void CurvePoints(int cur_num,int totalpoints,vector <int> curve_length,vector<Point> curve_set,int flag);

extern void drawpoints(vector<Point> Points);
extern void drawfeaturepoints(vector<Point> Points);
extern vector<Point> CurveSampling(int cur_num,int totalpoints,vector <int> &curve_length,
								   vector<Point> curve_set,int flag,vector <int> &start,vector <int>&end);
#endif