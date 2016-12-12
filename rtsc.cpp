#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include "GL/glui.h"
#include "TriMesh.h"
#include "TriMesh_algo.h"
#include "XForm.h"
#include "GLCamera.h"
#include "timestamp.h"
#include "windows.h"
#include "Commdlg.h"
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include "Wire.h"
#include "math.h"
#include <iostream>

using std::min;
using std::max;
using std::swap;
#define BIG 10000000000000
#define SMALL -1e20
#define LITTLE 1e-6

#define N_Models 800
//#define NEEDSAMPLING 1
#define TWOPI 6.2831852
#define LEFT -1.0
#define RIGHT 1.0
#define BOTTOM -1.0
#define TOP 1.0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define Thresh 3

//State of Mouse pressing 
#define eps 1e-6

//icosahedron lonely variables
#define V_N0 12 //	6	//level 0  
#define V_N1 42//	21 	//level 1  
#define V_N2 162  //81	//level 2   
#define V_N3  642//	321//level 3
#define V_N4  2562  //1281	//level 4
#define V_N5 10242	//level 5
#define V_N6 40002 	//level 5
int VNUM=0;
point O;

/*added by libo,20080829,begin*/
int GW=400, GH=400;	//window size
//int GW=128,GH=128;

int  whichMethod, whichUpright, whichColor, whichBruteLevel, whichShading, bruteID=0, Metro=1;
float T_diff_eps = 1.5e-2; 
point pointNol, alignAxes[3];

vector<point> flock; 
double tmpArea[100000]={0};
float pixelColorBuf[400][400][3]={0}; 
static float AR;
float scale_ratio=1.0;//1.5; //control the ratio of projection,default [-1.5,1.5]
float  xUser, yUser, zUser;
GLUI *glui;

TriMesh *themesh, *Cmesh, *icoMesh0, *icoMesh1, *icoMesh2, *icoMesh3, *icoMesh4, *icoMesh5, *icoMesh6;
double minS=BIG, maxS=SMALL; 
//Mesh modelMesh;
GLCamera camera;
xform xf;
int LIGHTING= FALSE;
int main_win,second_win,third_win;
void disp2();
void disp3();
void need_redraw();
//void drawfeaturepoints(vector<Point> Points);
void mymouse(int button, int state, int x, int y);
void reshape1(int w, int h);
void setEnv(); 
void KeyFunc(unsigned char key, int x, int y);
int  modelfile();
void quicksort(float arr[],int beg,int end);
void quicksort_double(double arr[],int beg,int end);
void drawModel();
void drawIcos(TriMesh *a);
float calVisibleSaliency();
double calVisibleEntropy();
void mesh_cpca_sign();
int  ReadMesh(char* szFile);
void CPCA_Normalization();
void drawSymmetryPlane(vector<point> flock);


// Clear the screen and reset OpenGL modes to something sane
void cls()
{
	glDisable(GL_DITHER);
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_NORMALIZE);
	glDisable(GL_LIGHTING);
	glDisable(GL_COLOR_MATERIAL);
	glClearColor(1,1,1,0);
	glClearDepth(1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

// Draw the scene
void disp1()
{
	timestamp t1 = now();
	// Set up camera and clear the screen
	GLUI_Master.auto_set_viewport();
	camera.setupGL(O, 1);
	cls(); // Clear the screen and reset OpenGL modes to something sane
	point v;
   
	// Draw the mesh using the correct drawing style
	glPushMatrix();
	glMultMatrixd((double *)xf);

	 glPopMatrix();
	 glutSwapBuffers();
	 fflush(stdout);	
	 if (camera.autospin(xf))
		need_redraw();	
}

// Signal a redraw
void need_redraw()
{	
	glutPostRedisplay();
}

 	
float axis1[3]={1,0,0};  //rotate about X axis for phai
float axis2[3]={0,1,0};  //rotate about Y axis for theta
float axis3[3]={0,0,1};  //rotate about Z axis for theta
// Set the view to look at the middle of the mesh, from reasonably far away
int pose=0;

void resetView()
{
	xf=xform::trans(0,0,0);
	//xf = xf*xform::trans(0, 0, -10.0f * themesh->bsphere.r);
	camera.stopspin();
}

void recoverPoseXYZ(float x, float y, float z)
{
	resetView();
	float r=sqrt(x*x+y*y+z*z);	
	float angle=acos(z/r);
	xf=xf*xform::rot(angle,y/r,-1*x/r,0);
	glutSetWindow(second_win);
	disp2();
}

void readIcos()
{
	int i; vec2 sc;
	Color c = Color::black();

	//read icosahedron with depth level 0,1,2,3,and 4 (12, 42, 162, 642, and 2562 vertices)
	icoMesh0 = TriMesh::read("..\\icos\\icos0.off");
	icoMesh1 = TriMesh::read("..\\icos\\icos1.off");
	icoMesh2 = TriMesh::read("..\\icos\\icos2.off");
	icoMesh3 = TriMesh::read("..\\icos\\icos3.off");
	icoMesh4 = TriMesh::read("..\\icos\\icos4.off");
	icoMesh5 = TriMesh::read("..\\icos\\icos5.off");
	icoMesh6 = TriMesh::read("..\\icos\\icos6.off");	
	//create neighbor data structure and others
	icoMesh0->need_neighbors();
	icoMesh0->need_faces();
	icoMesh0->need_normals();
	for (i=0;i<(int)icoMesh0->vertices.size();i++)
		icoMesh0->colors.push_back(c);

	icoMesh1->need_neighbors();
	icoMesh1->need_faces();
	icoMesh1->need_normals();
	for (i=0;i<(int)icoMesh1->vertices.size();i++)
		icoMesh1->colors.push_back(c);

	icoMesh2->need_neighbors();
	icoMesh2->need_faces();
	icoMesh2->need_normals();
	for (i=0;i<(int)icoMesh2->vertices.size();i++)
		icoMesh2->colors.push_back(c);

	icoMesh3->need_neighbors();
	icoMesh3->need_faces();
	icoMesh3->need_normals();
	for (i=0;i<(int)icoMesh3->vertices.size();i++)
		icoMesh3->colors.push_back(c);

	icoMesh4->need_neighbors();
	icoMesh4->need_faces();
	icoMesh4->need_normals();
	for (i=0;i<(int)icoMesh4->vertices.size();i++)
		icoMesh4->colors.push_back(c);

	icoMesh5->need_neighbors();
	icoMesh5->need_faces();
	icoMesh5->need_normals();
	for (i=0;i<(int)icoMesh5->vertices.size();i++)
		icoMesh5->colors.push_back(c);

	icoMesh6->need_neighbors();
	icoMesh6->need_faces();
	icoMesh6->need_normals();
	for (i=0;i<(int)icoMesh6->vertices.size();i++)
		icoMesh6->colors.push_back(c);	
}

 void quicksort_double(double arr[],int beg,int end)
{
	if (end  >= beg + 1) 
  	{
  		double piv = arr[beg];
		int k = beg + 1, r = end;
    
		while (k < r) 
    		{
      			if (arr[k] > piv) 
        			k++;
      			else 
        			swap(arr[k], arr[r--]);
    		}
		if (arr[k] > piv){
		
			swap(arr[k],arr[beg]);
			
			quicksort_double(arr, beg, k);
			quicksort_double(arr, r, end);			
		}else {
			if (end - beg == 1)
  				return;
  				
			swap(arr[--k],arr[beg]);
			quicksort_double(arr, beg, k);
			quicksort_double(arr, r,   end);			
		}
  	}
}

void quicksort(float arr[],int beg,int end)
{
	if (end  >= beg + 1) 
  	{
  		float piv = arr[beg];
		int k = beg + 1, r = end;
    
		while (k < r) 
    		{
      			if (arr[k] < piv) 
        			k++;
      			else 
        			swap(arr[k], arr[r--]);
    		}
		if (arr[k] < piv){
		
			swap(arr[k],arr[beg]);
			
			quicksort(arr, beg, k);
			quicksort(arr, r, end);			
		}else {
			if (end - beg == 1)
  				return;
  				
			swap(arr[--k],arr[beg]);
			quicksort(arr, beg, k);
			quicksort(arr, r,   end);			
		}
  	}
}

void idQuickSort(int id[], double arr[], int beg, int end) {
	if (end <= beg) {
		return;
	}
	int i = beg;
	int j = beg;
	while (j < end) {
		if (arr[id[j]] < arr[id[end]]) {
			if (i != j) {
				swap(id[i], id[j]);
			}
			++i;
		}
		++j;
	}
	swap(id[i], id[end]);
	idQuickSort(id, arr, beg, i - 1);
	idQuickSort(id, arr, i + 1, end);
}

// Handle mouse button and motion events
static unsigned buttonstate = 0;

void mousemotionfunc(int x, int y)
{
	static const Mouse::button physical_to_logical_map[] = {
		Mouse::NONE, Mouse::ROTATE, Mouse::MOVEXY, Mouse::MOVEZ,
		Mouse::MOVEZ, Mouse::MOVEXY, Mouse::MOVEXY, Mouse::MOVEXY,
	};
	Mouse::button b = Mouse::NONE;
	if (buttonstate & (1 << 3))
		b = Mouse::WHEELUP;
	else if (buttonstate & (1 << 4))
		b = Mouse::WHEELDOWN;
	else
		b = physical_to_logical_map[buttonstate & 7];
	//reset here, added by libo
	if (b == Mouse::ROTATE)
	{	
	}
     
	camera.mouse(x, y, b,
		     O, 1,
		     xf);

      	 /*
	xf.rot(90/180,0,0,1);
	vec axis(0,0,1);
	rot(themesh,30/180,axis);
       */
	need_redraw();
}

void mousebuttonfunc(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN)
		buttonstate |= (1 << button);
	else
		buttonstate &= ~(1 << button);

	mousemotionfunc(x, y);	
}


// Keyboard callback
void keyboardfunc(unsigned char key, int x, int y)
{
	switch (key) {
		case 'Q':
		case 'q':
			exit(0);		
	}
	need_redraw(); 
	GLUI_Master.sync_live_all();
}

void skeyboardfunc(int key, int x, int y)
{
	need_redraw(); 
	GLUI_Master.sync_live_all();
}

void reshape(int x, int y)
{
	GLUI_Master.auto_set_viewport();
	cls();
	glutSwapBuffers();
	need_redraw();
}


void usage(const char *myname)
{
	fprintf(stderr, "Usage: %s [-options] infile\n", myname);
	exit(1);
}
double calSymmetry()
{
	int i, j, k;
	double area=0,m=0;

	//addition
	float centerX=0,centerY=0, areaY=0, areaX=0;
	int x,y;

	int View[4];
	glutSetWindow(second_win);
	glGetIntegerv(GL_VIEWPORT, View);
	int width = View[2], height = View[3];
	centerX=width/2;
	centerY=height/2;

	for (i = 0; i < width; i++)
	for (j = 0; j < height; j++)
	  for (k = 0; k < 3; k++)
		  pixelColorBuf[i][j][k] = 0;
	glReadBuffer(GL_FRONT);
	glReadPixels(View[0], View[1], width, height, GL_RGB, GL_FLOAT, pixelColorBuf);
	//maximize symmetry along y axis
	//check left to right
	for(i=0;i<width/2;i++)
	{
		for (j=0;j<height;j++)
		{
			if (pixelColorBuf[i][j][2]==0.0)
			{
				x=width-i;
				if(x!=i)
				if(pixelColorBuf[x][j][2]==0.0)
					areaY+=log(centerX-i);
			}
		}
	}
	//maximize symmetry along x axis
	//check up to bottom
	for(i=0;i<width;i++)
	{
		for (j=0;j<height/2;j++)
		{
			if (pixelColorBuf[i][j][2]==0.0)
			{
				 y=height-j;
				 if(y!=j)
					 if(pixelColorBuf[x][j][2]==0.0)
						areaX+=log(centerY-j);
			}
		}
	}

	area = areaX;//areaY; //make the object view as tall as possible but as slim as possible by playing around with symmetry property
	return area;
}

double calArea(int v)
{
	int i, j, k;
	double area=0,m=0;

	int View[4];
	glutSetWindow(second_win);
	glGetIntegerv(GL_VIEWPORT, View);
	int width = View[2], height = View[3];

	for (i = 0; i < width; i++)
		for (j = 0; j < height; j++)
			for (k = 0; k < 3; k++)
				pixelColorBuf[i][j][k] = 0;
	glReadBuffer(GL_FRONT);
	glReadPixels(View[0], View[1], width, height, GL_RGB, GL_FLOAT, pixelColorBuf);
	for (i = 0; i < width; i++){
		for (j = 0; j < height; j++)  
			if (fabs(pixelColorBuf[i][j][2])<eps) 
				area=area+1;
	}
	//area=area/(GW*GH);
	if(bruteID==1)
		tmpArea[v]=area;
	return area;
}
	
void ComputeViewFeatures()
{
	int i;
	double s;
	point p;
	for(i=0;i<(int)Cmesh->vertices.size();i++)
	{
		p = Cmesh->vertices[i];						
		recoverPoseXYZ(p[0],p[1],p[2]);
		switch(whichUpright)
		{
			case 0: 
				s=calArea(i); 
				break;
			case 1:
				s=calVisibleEntropy();					
				break;
			case 3: 
				s=calSymmetry(); 
				break;
		}
		tmpArea[i]=s;
		if (minS>s)
			minS=s;
		if (maxS<s)
			maxS=s;
	}
}

void computeColorDist()
{
	//gray and rgb variables
	double range = maxS-minS, tmp, tmpRGB, maxR=1, maxG=2, maxB =3;
	int i;
	//hsv variables
	float hue, tmpHue, v=1, s=1, X, r, g, b;
	int index;
	
	for(i=0;i<VNUM;i++)
	{
		tmp=(tmpArea[i]-minS)/range;
		tmpRGB=tmp*3;
		for(int k=0;k<3;k++) Cmesh->colors[i][k]=0;
		//map areas to colors
		if(whichColor==0)//grayscale colors
			for(int j=0;j<3;j++) Cmesh->colors[i][j]=tmp;
		else if(whichColor==1)
		{//rgb colors
			if(tmpRGB<maxR)Cmesh->colors[i][0]=tmpRGB;
			else if(tmpRGB>=maxR && tmpRGB<maxG)Cmesh->colors[i][1]=tmpRGB-1;
			else Cmesh->colors[i][2]=tmpRGB-2;
		}
		else
		{//map to HSV for natural color then remap back to rgb
			hue=(tmpArea[i]-minS)/range*360; //map to 0-360
			tmpHue=(hue*330/360)-60; //map to (-60 to 270)
			if(tmpHue<0) tmpHue+=60;
			//remap into rgb space
			tmpHue/=60;			// sector 0 to 5
			index = floor(tmpHue);
			X=v*s*(1-fabs(fmod(tmpHue,2)-1));
			switch( index ) {
				case 0:
					r = v*s; g = X; b = 0;break;
				case 1:
					r = X;g = v*s;b = 0;break;
				case 2:
					r = 0;g = v*s;b = X;break;
				case 3:
					r = 0;g = X;b = v*s;break;
				case 4:
					r = X;g = 0;b = v*s;break;
				default:		//case 5:
					r = v*s;g = 0;b = X;break;
			}
			Cmesh->colors[i][0]=r;
			Cmesh->colors[i][1]=g;
			Cmesh->colors[i][2]=b;
		}
	}
}

float uaxis[3]={0};
float vaxis[3]={0};
float waxis[3]={0};
float mean[3]={0}, cov[3][3]={0};
#define SIGN(a, b) ( (b) < 0 ? -fabs(a) : fabs(a) )

/**  Error handler  **************************************************/
void erhand1(char err_msg[])
/* Error handler */
{
    fprintf(stderr,"Run-time error:\n");
    fprintf(stderr,"%s\n", err_msg);
    fprintf(stderr,"Exiting to system.\n");
    exit(1);
}
/**  Allocation of vector storage  ***********************************/

double* vector1(int n)
/* Allocates a float vector with range [1..n]. */
{

    double *v;

    v = (double *) malloc ((unsigned) n*sizeof(double));
    if (!v) 
		erhand1("Allocation failure in vector().");
	return v-1;

}
/**  Allocation of float matrix storage  *****************************/

double **matrix(int n,int m)
/* Allocate a float matrix with range [1..n][1..m]. */
{
    int i;
    double **mat;

    /* Allocate pointers to rows. */
    mat = (double **) malloc((unsigned) (n)*sizeof(double*));
    if (!mat) erhand1("Allocation failure 1 in matrix().");
    mat -= 1;

    /* Allocate rows and set pointers to them. */
    for (i = 1; i <= n; i++)
        {
        mat[i] = (double *) malloc((unsigned) (m)*sizeof(double));
        if (!mat[i]) erhand1("Allocation failure 2 in matrix().");
        mat[i] -= 1;
        }

     /* Return pointer to array of pointers to rows. */
     return mat;

}
/**  Reduce a real, symmetric matrix to a symmetric, tridiag. matrix. */

void tred2(double **a,int  n,double *d,double *e)
/* float **a, d[], e[]; */
/* Householder reduction of matrix a to tridiagonal form.
   Algorithm: Martin et al., Num. Math. 11, 181-195, 1968.
   Ref: Smith et al., Matrix Eigensystem Routines -- EISPACK Guide
        Springer-Verlag, 1976, pp. 489-494.
        W H Press et al., Numerical Recipes in C, Cambridge U P,
        1988, pp. 373-374.  */
{
int l, k, j, i;
double scale, hh, h, g, f;

for (i = n; i >= 2; i--)
    {
    l = i - 1;
    h = scale = 0.0;
    if (l > 1)
       {
       for (k = 1; k <= l; k++)
           scale += fabs(a[i][k]);
       if (scale == 0.0)
          e[i] = a[i][l];
       else
          {
          for (k = 1; k <= l; k++)
              {
              a[i][k] /= scale;
              h += a[i][k] * a[i][k];
              }
          f = a[i][l];
          g = f>0 ? -sqrt(h) : sqrt(h);
          e[i] = scale * g;
          h -= f * g;
          a[i][l] = f - g;
          f = 0.0;
          for (j = 1; j <= l; j++)
              {
              a[j][i] = a[i][j]/h;
              g = 0.0;
              for (k = 1; k <= j; k++)
                  g += a[j][k] * a[i][k];
              for (k = j+1; k <= l; k++)
                  g += a[k][j] * a[i][k];
              e[j] = g / h;
              f += e[j] * a[i][j];
              }
          hh = f / (h + h);
          for (j = 1; j <= l; j++)
              {
              f = a[i][j];
              e[j] = g = e[j] - hh * f;
              for (k = 1; k <= j; k++)
                  a[j][k] -= (f * e[k] + g * a[i][k]);
              }
         }
    }
    else
        e[i] = a[i][l];
    d[i] = h;
    }
d[1] = 0.0;
e[1] = 0.0;
for (i = 1; i <= n; i++)
    {
    l = i - 1;
    if (d[i])
       {
       for (j = 1; j <= l; j++)
           {
           g = 0.0;
           for (k = 1; k <= l; k++)
               g += a[i][k] * a[k][j];
           for (k = 1; k <= l; k++)
               a[k][j] -= g * a[k][i];
           }
       }
       d[i] = a[i][i];
       a[i][i] = 1.0;
       for (j = 1; j <= l; j++)
           a[j][i] = a[i][j] = 0.0;
    }
}

/**  Tridiagonal QL algorithm -- Implicit  **********************/

void tqli(double d[], double e[], int n, double **z)
{
int m, l, iter, i, k;
double s, r, p, g, f, dd, c, b;
void erhand();

for (i = 2; i <= n; i++)
    e[i-1] = e[i];
e[n] = 0.0;
for (l = 1; l <= n; l++)
    {
    iter = 0;
    do
      {
      for (m = l; m <= n-1; m++)
          {
          dd = fabs(d[m]) + fabs(d[m+1]);
          if (fabs(e[m]) + dd == dd) break;
          }
          if (m != l)
             {
             if (iter++ == 30) erhand1("No convergence in TLQI.");
             g = (d[l+1] - d[l]) / (2.0 * e[l]);
             r = sqrt((g * g) + 1.0);
             g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
             s = c = 1.0;
             p = 0.0;
             for (i = m-1; i >= l; i--)
                 {
                 f = s * e[i];
                 b = c * e[i];
                 if (fabs(f) >= fabs(g))
                    {
                    c = g / f;
                    r = sqrt((c * c) + 1.0);
                    e[i+1] = f * r;
                    c *= (s = 1.0/r);
                    }
                 else
                    {
                    s = f / g;
                    r = sqrt((s * s) + 1.0);
                    e[i+1] = g * r;
                    s *= (c = 1.0/r);
                    }
                 g = d[i+1] - p;
                 r = (d[i] - g) * s + 2.0 * c * b;
                 p = s * r;
                 d[i+1] = g + p;
                 g = c * r - b;
                 for (k = 1; k <= n; k++)
                     {
                     f = z[k][i+1];
                     z[k][i+1] = s * z[k][i] + c * f;
                     z[k][i] = c * z[k][i] - s * f;
                     }
                 }
                 d[l] = d[l] - p;
                 e[l] = g;
                 e[m] = 0.0;
             }
          }  while (m != l);
      }
 }
void SVD() 
{
	int i,j;
	double *interm,*evals,*tevals,**symmat;
	interm = vector1(3);    /* Storage alloc. for 'intermediate' vector */
	evals = vector1(3);     /* Storage alloc. for vector of eigenvalues */
 		
	symmat = matrix(3, 3);  /* Allocation of correlation (etc.) matrix */
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			symmat[i+1][j+1]=cov[i][j];

	tred2(symmat, 3, evals, interm);  /* Triangular decomposition */		
	tqli(evals, interm, 3, symmat);   /* Reduction of sym. trid. matrix */
	
	tevals = vector1(3);     /* Tempory storage alloc. for vector of eigenvalues */
	for(i=1;i<=3;i++)
		tevals[i]=evals[i];
	
	quicksort_double(tevals,1,3);
/*	printf("\nEigenvalues:\n");
	for (int j = 1; j <= 3; j++) 
	{
		printf("%18.6f\n", tevals[j]); 
	}
*/	
	int index[4]={0};
	for(i=1;i<=3;i++)
		for(j=1;j<=3;j++)		
			if(evals[j]==tevals[i])
			{
				index[i]=j;
				break;
			}
	for(i=0;i<3;i++)
		uaxis[i]=symmat[i+1][index[1]];
	for(i=0;i<3;i++)
		vaxis[i]=symmat[i+1][index[2]];	
	for(i=0;i<3;i++)
		waxis[i]=symmat[i+1][index[3]];
}
void meanFunc()
{
	double mx=0, my = 0, mz=0;
	for(int i=0;i<themesh->vertices.size();i++)
	{
	   mx=mx+themesh->vertices[i][0];
	   my=my+themesh->vertices[i][1]; 
	   mz=mz+themesh->vertices[i][2]; 
	}
	mean[0]=mx/themesh->vertices.size();
	mean[1]=my/themesh->vertices.size();
	mean[2]=mz/themesh->vertices.size();
}

void covariance()
{
    double x,y,z;
	int row=themesh->vertices.size();
	
	for(int i = 0;i<3;i++)
		for(int j= 0;j<3;j++)
			cov[i][j]=0;
		
	for(int i = 0;i<row;i++)
	{
	    x=themesh->vertices[i][0];
		y=themesh->vertices[i][1];
		z=themesh->vertices[i][2]; 	      
		cov[0][0]+=pow(x-mean[0],2);
		cov[0][1]+=(x-mean[0])*(y-mean[1]);
		cov[0][2]+=(x-mean[0])*(z-mean[2]);
		cov[1][1]+=pow(y-mean[1],2);
		cov[1][2]+=(y-mean[1])*(z-mean[2]);
		cov[2][2]+=pow(z-mean[2],2);cov[1][2]+=(y-mean[1])*(z-mean[2]);
	}

	
	cov[0][0]=cov[0][0]/row;
	cov[0][1]=cov[0][1]/row;
	cov[0][2]=cov[0][2]/row;
	cov[1][1]=cov[1][1]/row;
	cov[1][2]=cov[1][2]/row;
	cov[2][2]=cov[2][2]/row;
	
	
	cov[1][0]=cov[0][1];
	cov[2][0]=cov[0][2];
	cov[2][1]=cov[1][2];
	
	cout<<endl<<"Covariance matix: "<< cov[0][0]<<" "<<cov[0][1]<<" "<<cov[0][2]<<endl;
	cout<<"                  "<< cov[1][0]<<" "<<cov[1][1]<<" "<<cov[1][2]<<endl;
	cout<<"                  "<< cov[2][0]<<" "<<cov[2][1]<<" "<<cov[2][2]<<endl;
}
void CPCAcovariance()
{    
	mesh_cpca_covariance(themesh,cov);

//-	cout<<endl<<"Covariance matix: "<< cov[0][0]<<" "<<cov[0][1]<<" "<<cov[0][2]<<endl;
//-	cout<<"                  "<< cov[1][0]<<" "<<cov[1][1]<<" "<<cov[1][2]<<endl;
//-	cout<<"                  "<< cov[2][0]<<" "<<cov[2][1]<<" "<<cov[2][2]<<endl;
}
void CPCAmeanFunc()
{
       point com = mesh_center_of_mass(themesh);
   	mean[0]=com[0];
	mean[1]=com[1];
	mean[2]=com[2];
}

void dataprojection()
{
   float x,y,z;
   CPCAmeanFunc(); 
   for(int i=0;i<themesh->vertices.size();i++)
   {   
	   x=themesh->vertices[i][0];
	   y=themesh->vertices[i][1];
	   z=themesh->vertices[i][2]; 
	  
	   themesh->vertices[i][0]=uaxis[0]*(x-mean[0])+uaxis[1]*(y-mean[1])+uaxis[2]*(z-mean[2]);
	   themesh->vertices[i][1]=vaxis[0]*(x-mean[0])+vaxis[1]*(y-mean[1])+vaxis[2]*(z-mean[2]);
	   themesh->vertices[i][2]=waxis[0]*(x-mean[0])+waxis[1]*(y-mean[1])+waxis[2]*(z-mean[2]);
   }
}

float J(float xa,float xb,float xc)
{
   return (xa*xa+xb*xb+xc*xc+xa*xb+xa*xc+xb*xc);
}
float L(float xa,float xb,float xc)
{
  return(pow(xa,4)/((xb-xa)*(xc-xa)));
}
float F(TriMesh *mesh,int index, int axis)
{
   const TriMesh::Face &f = mesh->faces[index];
   const vector<point> &p = mesh->vertices;
   float coord[3]={0};
   int i=0;
   for (int v = 0; v < 3; v++) 
   {
       point pm = p[f[v]];
	   coord[i]=pm[axis];
	   i++;
   }
   quicksort(coord,0,2);
   if(coord[0]>=0)
   	return (J(coord[0],coord[1],coord[2]));
   else if(coord[2]<0)
   	return (-J(coord[0],coord[1],coord[2]));
   else if(coord[0]<0 && coord[1]>=0 && coord[2]>=0)
   	return (J(coord[0],coord[1],coord[2])-2*L(coord[0],coord[1],coord[2]));
   else if(coord[2]>=0 && coord[1]<0 && coord[0]<0)
   	return (-J(coord[2],coord[1],coord[0])+2*L(coord[2],coord[1],coord[0])); 	   	
  }
	
// Compute cpca covariance of faces (integral area-weighted) in a mesh
void mesh_cpca_sign()
{
	float Sign[3]={0};
	//point m = mesh_center_of_mass(themesh);
		
	int n = themesh->faces.size();
	
	// Compute area of each face and the total area
	float *area;
	area = new float[n];	
	const vector<point> &p = themesh->vertices;	
	for (int i = 0; i < n; i++) {
		const TriMesh::Face &f = themesh->faces[i];
		area[i] = len(trinorm(p[f[0]], p[f[1]], p[f[2]]));
	}
	
	
	for (int i = 0; i < n; i++) 
		for (int j = 0; j< 3; j++) 
		{		
			Sign[j] += area[i] * F(themesh,i,j);
		}			

	for (int j = 0; j < 3; j++)
	  if (Sign[j]<0)
	  	Sign[j]=-1;
	  else
	  	Sign[j]=1;	
	for(int i=0;i<themesh->vertices.size();i++) 
	   for (int j = 0; j< 3; j++) 
	   	themesh->vertices[i][j]=themesh->vertices[i][j]*Sign[j];
	delete [] area;
}
void CPCA_Normalization()
{ 
	CPCAmeanFunc();
	CPCAcovariance();
	SVD();	
	dataprojection();
	mesh_cpca_sign();
}

void output(TriMesh* mesh, char* filename)
{
	FILE *fp;
	int i;
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("Cannot open file\n");
		exit(0);
	}
	char head[10]="OFF";
	fprintf(fp,"%s",head);
	fprintf(fp,"\n");
	fprintf(fp,"%d %d %d\n",mesh->vertices.size(),mesh->faces.size(),0);

	for(i=0;i<mesh->vertices.size();i++)
	{
		fprintf(fp,"%f %f %f\n",mesh->vertices[i][0],mesh->vertices[i][1],mesh->vertices[i][2]);       
	}

	for (i=0;i<(int)mesh->faces.size();i++)
	{
		fprintf(fp,"%d %d %d %d\n",3,mesh->faces[i][0],mesh->faces[i][1],mesh->faces[i][2]);         
	}

	fclose(fp);
}

void viewPointBrute()
{	
	timestamp t1 = now(); //initiate time counter
	TriMesh *resultmesh;
	
	whichUpright=1;
	LIGHTING=0;	
	themesh->need_tstrips();
	CPCA_Normalization();
	themesh->need_bsphere();
	float r=themesh->bsphere.r;
	vec c=themesh->bsphere.center;
	if (Metro)
		resultmesh=new TriMesh();
	for(int i=0;i<themesh->vertices.size();i++)
	{   
	   themesh->vertices[i]=(themesh->vertices[i]-c)/r;
	   if (Metro)	  
	   {
		   resultmesh->vertices.push_back(themesh->vertices[i]);
	   }
	}
	
	if (Metro)
	{
		for(int i=0;i<themesh->faces.size();i++)
			resultmesh->faces.push_back(themesh->faces[i]);
		output(themesh, "cpca.off");
	}

	ComputeViewFeatures();  //  1
	computeColorDist();
	resetView();

	const float P_diff_eps = 1e-6;
	int indices = 0;
	
	int maxc = 0;

	bool marked[3] = {false, false, false};
	int u, v, i;
//	#pragma omp parallel for default(none)  private (u, v, i) shared(VNUM, Cmesh, T_diff_eps, tmpArea)
	for (u = 0; u < VNUM-1; ++u) {
		vec pu = Cmesh->vertices[u];
		for (v =u+1; v < VNUM; ++v) {
			if (fabs(tmpArea[u]- tmpArea[v]) > T_diff_eps*min(tmpArea[u],tmpArea[v])) {
				continue;
			}
			vec pv = Cmesh->vertices[v];
			vec T1=pu-pv;
			T1=normalize(T1);

			int counter = 2;
//			#pragma omp parallel for default(none)  private (i, u, v) shared(VNUM, Cmesh, tmpArea, T_diff_eps, T1, counter)
			for (i = 0; i < VNUM-1; ++i) {
				if (i == u || i == v) {
					continue;
				}
				vec pi = Cmesh->vertices[i];

				for (int j = i+1; j < VNUM; ++j) {
					if (j == u || j == v) {
						continue;
					}
					if (fabs(tmpArea[i] - tmpArea[j]) > T_diff_eps*min(tmpArea[i],tmpArea[j])) {
						continue;
					}
					vec pj = Cmesh->vertices[j];
					vec m=pi+pj;
					m*=0.5;
					
					vec T2 =pj-pi;
					T2=normalize(T2);
					
					vec CT=T1 CROSS T2;
					float DT=T1 DOT T2;

					if (len(CT)>P_diff_eps && fabs(DT)!=0){
						continue;
					}

					if (fabs(T1 DOT m) > P_diff_eps) {
						continue;
					}	
					
					counter=counter+2;
					break;
				}
			}
			
			//printf("counter = %d\n", counter);
			maxc = max(maxc, counter);
			float a=T1[0], b=T1[1], c=T1[2];			
			if (counter>=VNUM-int(pow(2.0, double(whichBruteLevel+2)))) {
				printf("Counter:%d:%f * x + %f * y +  %f * z = 0\n", counter,a, b, c);
			
				point p1;
				point p2;
				point p3;
				point p4;
				//if (a*b*c!=0)  //a, b, c !=0
				float L=0.98;
				if (fabs(a)>eps && fabs(b)>eps && fabs(c)>eps)  //a, b, c !=0
				{
					p1 = point(L, L, (a * (0 - L) + b * (0 - L)) / c);
					p2 = point(-L, L, (a * (0 + L) + b * (0 - L)) / c);
					p3 = point(L, -L, (a * (0 - L) + b * (0 + L)) / c);
					p4 = point(-L, -L, (a * (0 + L) + b * (0 + L)) / c);
				}
				else if (b*c!=0)  //only a=0
				{
				   	p1 = point(L, L, -b/c);
				   	p2 = point(L, -L, b/c);				   	
					p3 = point(-L, L, -b/c);			
					p4 = point(-L, -L,  b/c);
				}		
				else if (a*c!=0) //only b=0
				{
				   	p1 = point(L, L, -a/c);
				   	p2 = point(L, -L, a/c);				   	
					p3 = point(-L, L, -a/c);			
					p4 = point(L, -L, a/c);
				}					
				else if (a*b!=0) //only c=0
				{
				   	p1 = point(L, -a/b, L);
					p2 = point(-L, a/b, L);
				   	p3 = point(L, -a/b, -L);
					p4 = point(-L, a/b, -L);					
				}				
				else if (fabs(a)>eps)  //only a!=0
				{
				   	p1 = point(0, L, L);
				   	p2 = point(0, L, -L);				   	
					p3 = point(0, -L, L);			
					p4 = point(0, -L, -L);
				}
				else if (fabs(b)>eps)  //only b!=0
				{
				   	p1 = point(L, 0, L);
				   	p2 = point(L, 0, -L);				   	
					p3 = point(-L, 0, L);			
					p4 = point(-L, 0, -L);
				}	
				else if (fabs(c)>eps)  //only c!=0
				{
				   	p1 = point(L, L, 0);
				   	p2 = point(L, -L, 0);				   	
					p3 = point(-L, L, 0);			
					p4 = point(-L, -L, 0);
				}				
				flock.push_back(p1);
				flock.push_back(p2);
				flock.push_back(p3);
				flock.push_back(p4);
				if (Metro)
				{
					if (fabs(a) ==1) {
						if (false == marked[0]) {
							marked[0] = true;
							for (i = 0; i < resultmesh->vertices.size(); ++i) {
								resultmesh->vertices[i][0] = -resultmesh->vertices[i][0];
							}
							output(resultmesh, "reflection_x.off");
							//change back to keep the original
							for (int i = 0; i < resultmesh->vertices.size(); ++i) {
								resultmesh->vertices[i][0] = -resultmesh->vertices[i][0];
							}
						}
					} else if (fabs(b) ==1) {
						if (false == marked[1]) {
							marked[1] = true;
							for (i = 0; i < resultmesh->vertices.size(); ++i) {
								resultmesh->vertices[i][1] = -resultmesh->vertices[i][1];
							}
							output(resultmesh, "reflection_y.off");
							//change back to keep the original
							for (i = 0; i < resultmesh->vertices.size(); ++i) {
								resultmesh->vertices[i][1] = -resultmesh->vertices[i][1];
							}	
							//output(resultmesh, "original.off");
						}
					} else if (fabs(c) ==1) {
						if (false == marked[2]) {
							marked[2] = true;
							for (i = 0; i < resultmesh->vertices.size(); ++i) {
								resultmesh->vertices[i][2] = -resultmesh->vertices[i][2];
							}
							output(resultmesh, "reflection_z.off");
							//change back to keep the original							
							for (int i = 0; i < resultmesh->vertices.size(); ++i) {
								resultmesh->vertices[i][2] = -resultmesh->vertices[i][2];
							}							
						}
					}
				}
			}
		}
	}
	printf("VNUM = %d, maxc = %d\n", VNUM, maxc);
	whichUpright=0;
	LIGHTING=1;
	resetView();
	glutSetWindow(second_win);
	disp2();
	glutSetWindow(third_win);	
	disp3();
	cout<<endl<<"Time: ...."<<now()-t1<<endl;    
}

void visualize()
{
	flock.clear();
	glutPostRedisplay();		
	switch (whichBruteLevel)
	{
		case 0: 
			VNUM=V_N0;
			Cmesh=icoMesh0;
			break;
		case 1: 
			VNUM=V_N1;
			Cmesh=icoMesh1;
			break;
		case 2: 
			VNUM=V_N2;
			Cmesh=icoMesh2;
			break;
		case 3: 
			VNUM=V_N3;
			Cmesh=icoMesh3;
			break;
		case 4: 
			VNUM=V_N4;
			Cmesh=icoMesh4;			
			break;
		case 5: 
			VNUM=V_N5;
			Cmesh=icoMesh5;				
			break;
		case 6: 
			VNUM=V_N6;
			Cmesh=icoMesh6;				
			break;	
	}
			
	for(int i=0;i<VNUM;i++) 
		tmpArea[i]=0;			

	modelfile();
	viewPointBrute();
}

char ModelName[N_Models][100];
char PureModelName[N_Models][50];
void ReadModelName()
{
	int i;

	 char modelfilename[256]="..\\NIST\\NistModelList.txt";
	   
	ifstream fin(modelfilename);
	for(i=0;i<N_Models;i++)
	{
		char name[50];
		fin.getline(name,128); 
		strcpy(ModelName[i],name);
		strcpy(PureModelName[i],name);
	}
}

void rotateToHere()
{
	float mag=sqrt(pow(xUser,2)+pow(yUser,2)+pow(zUser,2));
	if (mag==0)	xUser=LITTLE;
	xUser/=mag;
	yUser/=mag;
	zUser/=mag;
	glui->sync_live();
	LIGHTING=FALSE;
	recoverPoseXYZ(xUser, yUser, zUser);
	glutSetWindow(second_win);
	glui->sync_live();
//-	LIGHTING=TRUE;
	disp2();
}


int main(int argc, char *argv[])
{  
	//initialize the GLUT
	glutInitWindowSize(GW+210,GH+400);            
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInit(&argc, argv);

	if (argc < 2)
		usage(argv[0]);

	// Skip over any parameter beginning with a '-'
	int i = 1;
	while (i < argc-1 && argv[i][0] == '-')
		i++;
	const char *filename = argv[i];
	themesh = TriMesh::read(filename);
	readIcos();
	if (!themesh)
		usage(argv[0]);
 
	themesh->need_tstrips();
	themesh->need_bsphere();
	themesh->need_normals();

	float r=themesh->bsphere.r;
	vec c=themesh->bsphere.center;
	for(int i=0;i<themesh->vertices.size();i++)
	{   
	   themesh->vertices[i]=(themesh->vertices[i]-c)/r;
	}
	Cmesh = new TriMesh();
	
   	char windowname[255];
	sprintf(windowname, "RTSC - %s", filename);
	glutInitWindowPosition(200,100);
	main_win = glutCreateWindow(windowname);
	

	glutDisplayFunc(disp1);
	GLUI_Master.set_glutMouseFunc(mousebuttonfunc);
	glutMotionFunc(mousemotionfunc);
	GLUI_Master.set_glutKeyboardFunc(keyboardfunc);
	GLUI_Master.set_glutSpecialFunc(skeyboardfunc);
	GLUI_Master.set_glutReshapeFunc(reshape);

	//GLUI *glui = GLUI_Master.create_glui_subwindow(main_win, GLUI_SUBWINDOW_BOTTOM);
	glui = GLUI_Master.create_glui_subwindow(main_win, GLUI_SUBWINDOW_BOTTOM);
	glui->set_main_gfx_window(main_win);
	
	//initialize radio button state defaults
	whichShading=0;
	whichBruteLevel=1;
	whichColor=2; whichMethod=1,whichUpright=1; 
	
	GLUI_Panel *h = glui->add_panel("3D Symmetry Detection",1);
	//glui->add_button_to_panel(h ,"Load",0,(GLUI_Update_CB)modelfile);

	glui->add_separator_to_panel(h);
	glui->add_statictext_to_panel(h,"View Feature Selection:");
	GLUI_RadioGroup *up = glui->add_radiogroup_to_panel(h,&whichUpright);
	glui->add_radiobutton_to_group(up, "Area");  //MultiLevel
	glui->add_radiobutton_to_group(up, "Entropy");
	glui->add_separator_to_panel(h);
	//glui->add_edittext_to_panel(h, "Threshold:", GLUI_EDITTEXT_FLOAT, &T_diff_eps);	
	GLUI_Spinner *Tdiff= glui->add_spinner_to_panel(h,"Threshold",3,&T_diff_eps);
	Tdiff->set_float_limits(0,1,1);
	glui->add_separator_to_panel(h);
	
	glui->add_column_to_panel(h, 1);
	glui->add_separator_to_panel(h);
	glui->add_statictext_to_panel(h,"Level~Vertices");
	GLUI_RadioGroup *t = glui->add_radiogroup_to_panel(h,&whichBruteLevel);
		glui->add_radiobutton_to_group(t, "L0 ~ 12");
		glui->add_radiobutton_to_group(t, "L1 ~ 42");
		glui->add_radiobutton_to_group(t, "L2 ~ 162");
		glui->add_radiobutton_to_group(t, "L3 ~ 642");
		glui->add_radiobutton_to_group(t, "L4 ~ 2562");
		glui->add_radiobutton_to_group(t, "L5 ~ 10242");
		glui->add_radiobutton_to_group(t, "L6 ~ 40002");	
       glui->add_separator_to_panel(h);
	   
	glui->add_column_to_panel(h, 2);	   	
       glui->add_separator_to_panel(h);	
	glui->add_statictext_to_panel(h,"Shading");
	GLUI_RadioGroup *u = glui->add_radiogroup_to_panel(h,&whichShading);
		glui->add_radiobutton_to_group(u, "Smooth");
		glui->add_radiobutton_to_group(u, "Flat");
			
       glui->add_separator_to_panel(h);	
	glui->add_statictext_to_panel(h,"Metro Evaluation");
	GLUI_RadioGroup *m = glui->add_radiogroup_to_panel(h,&Metro);
		glui->add_radiobutton_to_group(m, "No");
		glui->add_radiobutton_to_group(m, "Yes");	
	glui->add_separator_to_panel(h);
	
	glui->add_column_to_panel(h, 3);			
	glui->add_separator_to_panel(h);	
	glui->add_button_to_panel(h ,"Detection",0,(GLUI_Update_CB)visualize); 

	glui->add_separator_to_panel(h);	
       glui->add_button_to_panel(h , "Exit", 0, exit);
	glui->add_separator_to_panel(h);	   
	resetView();

	//Go through command-line arguments and do what they say.
	for (int i = 1; i < argc-1; i++) 
	{
		if (argv[i][0] != '-')
			break;
		for (unsigned int j = 1; j < strlen(argv[i]); j++)
			keyboardfunc(argv[i][j], 0, 0);
	}
	
	//create the second  window for model display
		glutInitWindowSize(GW,GH);
		glutInitWindowPosition(725,100);
		second_win=glutCreateWindow("Sketch");
//		printf("second_win=%d\n",second_win);
		glViewport(0,0,GW,GH );
		//set the clearcolor and the callback
		glClearColor(0.0,0.0,0.0,1.0);
		glutDisplayFunc(disp2);
		glutReshapeFunc(reshape1);
		//glutMouseFunc(mymouse);
		GLUI_Master.set_glutMouseFunc(mousebuttonfunc);
		glutMotionFunc(mousemotionfunc);
		GLUI_Master.set_glutKeyboardFunc(keyboardfunc);
		GLUI_Master.set_glutSpecialFunc(skeyboardfunc);
		GLUI_Master.set_glutReshapeFunc(reshape);		
		glutKeyboardFunc(KeyFunc);
		//create the third window for visualization of cost function distribution
		glutInitWindowPosition(725,435);
		third_win=glutCreateWindow("Distribution");
//		printf("third_win=%d\n",third_win);
		glViewport(0,0,GW,GH );
		//set the clearcolor and the callback
		glClearColor(0.0,0.0,0.0,1.0);
		glutDisplayFunc(disp3);
		GLUI_Master.set_glutMouseFunc(mousebuttonfunc);
		glutMotionFunc(mousemotionfunc);
		glutReshapeFunc(reshape1);
		//glutMouseFunc(mymouse);
		glutKeyboardFunc(KeyFunc);
    glutMainLoop();	
 }
void drawModel()
{
	int i, j, id, num1, num2, index;
	point nn1,nn2,nn3;
	float tmpS, X, v=1,s=1, red,green,blue;
	glEnable(GL_DEPTH_TEST);
	point dummy, dummyN;	
	if (LIGHTING){
		glEnable(GL_LIGHTING); 	
		themesh->need_normals();
	}
	//glShadeModel(GL_FLAT);	
	glShadeModel(GL_SMOOTH);	
	xform nxf=norm_xf(xf);
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);
	int U=(int)sqrt((float)themesh->faces.size()); 
	for (i=0;i<(int)themesh->faces.size();i++)
	{
	    point v1 =xf* themesh->vertices[themesh->faces[i][0]];
		point v2 = xf*themesh->vertices[themesh->faces[i][1]];	
		point v3 = xf*themesh->vertices[themesh->faces[i][2]];
		if (LIGHTING)
		{
		    nn1=nxf*themesh->normals[themesh->faces[i][0]];  		
		    nn2=nxf*themesh->normals[themesh->faces[i][1]];
		    nn3 =nxf* themesh->normals[themesh->faces[i][2]];
		}
			
		glPointSize(1);	
		glBegin(GL_TRIANGLE_STRIP);
		if (whichUpright==1)   //entropy model
		{
			num1=i/U;
			num2=i-U*num1;
			glColor3f(num1/(float)(U-1),num2/(float)(U-1),0.0);
			glVertex3fv(scale_ratio*v1);			
			glVertex3fv(scale_ratio*v2);
			glVertex3fv(scale_ratio*v3);
			glEnd();
		}		
		else
		{
				glColor3f(0.0,0.0,0.0);		
				
				if (LIGHTING)
					glNormal3f(nn1[0], nn1[1], nn1[2]);				
				glVertex3fv(scale_ratio*v1);	

				if (LIGHTING)
					glNormal3f(nn2[0], nn2[1], nn2[2]);									
				glVertex3fv(scale_ratio*v2);

				if (LIGHTING)
					glNormal3f(nn3[0], nn3[1], nn3[2]);					
				glVertex3fv(scale_ratio*v3);
			
			}
		glEnd();
	} 
}
/*----------------------------------------------------------------------------
   Rendering routine
   --------------------------------------------------------------------------*/
void disp2() 
{  
   setEnv(); 
   glMatrixMode(GL_PROJECTION); 
   glLoadIdentity(); 
   glOrtho(-AR, AR, -1, 1, -100, 100);   
   glDisable(GL_LIGHTING);
   glClearColor (1.0, 1.0, 1.0, 0.0);
   if (themesh!=NULL)
   	drawModel();  
   if(!flock.empty())
   	drawSymmetryPlane(flock);   
   glutSwapBuffers();
   glutPostRedisplay();
}

 //Load A Model file
int  modelfile()
{
	int i,j;
	point centerMesh;
	OPENFILENAME OpenFileName;
	char szFile[MAX_PATH];

    //Create File Browser// 
	szFile[0] = 0;
	
    HWND hwnd= GetActiveWindow();              // owner window
    //printf("\nhwnd=%d\n",hwnd);
	OpenFileName.lStructSize = sizeof( OPENFILENAME );
	OpenFileName.hwndOwner = (HWND)hwnd;
	//OpenFileName.lpstrFilter =(LPCSTR)"\PLY(*.ply)\0*.ply\0All Files(.*)\0*.*\0";
	OpenFileName.lpstrFilter =(LPCSTR)"OFF(*.off)\0*.off\0All Files(.*)\0*.*\0";
	OpenFileName.lpstrCustomFilter = NULL;
	OpenFileName.nMaxCustFilter = 0;
	OpenFileName.nFilterIndex = 0;
	OpenFileName.lpstrFile = (LPSTR)szFile;
	OpenFileName.nMaxFile = sizeof( szFile );
	OpenFileName.lpstrFileTitle = NULL;
	OpenFileName.nMaxFileTitle = 0;
	OpenFileName.lpstrInitialDir = NULL;
	OpenFileName.lpstrTitle = (LPCSTR)"Open a file";
	OpenFileName.nFileOffset = 0;
	OpenFileName.nFileExtension = 0;
	OpenFileName.lpstrDefExt = NULL;
	OpenFileName.lCustData = 0;
	OpenFileName.lpfnHook = NULL;
	OpenFileName.lpTemplateName = NULL;
	OpenFileName.Flags = OFN_EXPLORER;

    //Load model 
	if (GetOpenFileName(&OpenFileName)==TRUE) 
	{
		delete themesh;
		themesh = TriMesh::read(szFile);
		if (themesh!=0)
		{		
			return TRUE;
		}
		else 
			return FALSE;
	}
	else 
	{   
		printf("Loading file  Failed!");
		return FALSE;
	}	
}
int  ReadMesh(char* szFile)
{
	int i,j;
    //Load model 
	{
		//printf("\n%s\n",szFile);	
		need_redraw();
		themesh = TriMesh::read(szFile);
		if (themesh!=0)
		{			
//			printf("Successfully loaded the model!");		
			themesh->need_tstrips();
			themesh->need_bsphere();
			//themesh->need_bbox();
			themesh->need_normals();

			resetView();
			//glutSetWindow(main_win);
			//disp1();
			glutPostRedisplay();
			return TRUE;
		}
		else 
			return FALSE;
	}
}


/*---------------------------------------------------------------------------
   Mouse Press p callback routine
   --------------------------------------------------------------------------*/
void mymouse(int button, int state, int x, int y)
{
	;
}

/*---------------------------------------------------------------------------
  Resize the window
  ---------------------------------------------------------------------------*/
void reshape1(int w, int h)
{
  GW = w;
  GH = h;
  AR = (float)(w)/(float)(h);
  glViewport(0, 0, w, h);               /* Set Viewport */   
  glutPostRedisplay();
}

/*---------------------------------------------------------------------------
  Set Rendering Enviroment
   --------------------------------------------------------------------------*/
void setEnv()
{
   glEnable(GL_DEPTH_TEST); 
   // Just clean the screen
   glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT); 
   
   // lighting
  // GLfloat light_position[]={ 1.0,1.0,1.0,0.0 };
   GLfloat light_position[]={ 0.0,0.0,1.0,0.0 };

   //GLfloat light_position[]={ 0.0,0.0,10.0,0.0 };
   GLfloat white_light[]={ 1.0,1.0,1.0,1.0 };
   GLfloat lmodel_ambient[]={ 1.0,0.1,0.1,1.0 };
   glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
   glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);

   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
   glEnable(GL_LIGHTING);glEnable(GL_LIGHT0);
   

   //Setup the camera
   glMatrixMode(GL_MODELVIEW); 
   glLoadIdentity(); 
   gluLookAt(0,0,20,0,0,0,0,1,0); 
   
    // setup the perspective projectin
   glMatrixMode(GL_PROJECTION); 
   //glLoadIdentity(); 
   glOrtho(-AR, AR, -1, 1, -100, 100);
}
/*---------------------------------------------------------------------------
   Key events function 
   --------------------------------------------------------------------------*/
void KeyFunc(unsigned char key, int x, int y)
{
  char sketch=0;
  switch(key)
  {      
 case 'q':
 case 'Q':
	 exit(0);
	 break;
  }
  glutPostRedisplay();
}

void disp3() 
{
	glLineWidth(6.0);
	setEnv(); 
	glMatrixMode(GL_PROJECTION); 
	glLoadIdentity(); 
	glOrtho(-AR*1.25, AR*1.25, -1.25, 1.25, -100, 100);  
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity(); 
	glDisable(GL_LIGHTING);
	glClearColor (1.0, 1.0, 1.0, 0.0);
	if (Cmesh!=NULL)
		drawIcos(Cmesh);
	if(!flock.empty())
		drawSymmetryPlane(flock);
	glutSwapBuffers();
	glutPostRedisplay();
}
void drawIcos(TriMesh *mesh)
{
	int i,j; 
	point nn1, nn2, nn3;
	glEnable(GL_DEPTH_TEST);
	glDisable (GL_LIGHTING);
	if(whichShading) 
		glShadeModel(GL_FLAT);	
	else 
		glShadeModel(GL_SMOOTH);

	xform nxf=norm_xf(xf);
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);

	for (i=0;i<(int)mesh->faces.size();i++)
	{
		point v1 =xf* mesh->vertices[mesh->faces[i][0]]; 
		point v2 = xf*mesh->vertices[mesh->faces[i][1]];
		point v3 = xf*mesh->vertices[mesh->faces[i][2]];

			
		glPointSize(1);
		glBegin(GL_TRIANGLES);	
			/***************first vertex*********************/
			glColor3fv(mesh->colors[mesh->faces[i][0]]);
			if (LIGHTING)
			{
			    nn1=nxf*mesh->normals[mesh->faces[i][0]]; 
				glNormal3f(nn1[0], nn1[1], nn1[2]);
			}
			glVertex3fv(scale_ratio*v1);			
	
			/***************second vertex*********************/
			glColor3fv(mesh->colors[mesh->faces[i][1]]);
			if (LIGHTING)
			{
				nn2=nxf*mesh->normals[mesh->faces[i][1]];
				glNormal3f(nn2[0], nn2[1], nn2[2]);
			}
				
			glVertex3fv(scale_ratio*v2);
			
			/***************third vertex*********************/
			glColor3fv(mesh->colors[mesh->faces[i][2]]);
			if (LIGHTING)
			{
				nn3 =nxf* mesh->normals[mesh->faces[i][2]];	
				glNormal3f(nn3[0], nn3[1], nn3[2]);
			}
			glVertex3fv(scale_ratio*v3);		
		glEnd();
			
	} 
}

void drawSymmetryPlane(vector<point> flock)
{
	int i, j, k;
	point v; 
	glLineWidth(1.0);
	glDisable(GL_LIGHTING);  

	glColor3f(143/255.0, 0, 1.0);
	for(j=0;j<(int)flock.size();j=j+4)
	{
		glBegin(GL_QUAD_STRIP);
		   for(k=j;k<j+4;k++)
		   {
		   	v=flock[k];
			glVertex3fv(xf*v);                
		    } 			 
		glEnd();
	}
}

float pixelcolorbuf[400][400][3] ={0}; 
double calVisibleEntropy()
{
  int i, j, k,index;
  double entropy=0.0;

  double *area,m=0,S;
  
  area=new double[themesh->faces.size()];
  for(i=0;i<themesh->faces.size();i++)
  	area[i]=0;
  
  int View[4];
  glutSetWindow(second_win);
  glGetIntegerv(GL_VIEWPORT, View);
  int width = View[2], height = View[3];

  for (i = 0; i < width; i++)
    for (j = 0; j < height; j++)
      for (k = 0; k < 3; k++)
		  pixelcolorbuf[i][j][k] = 0;
  
  glReadBuffer(GL_FRONT);
  glReadPixels(View[0], View[1], width, height, GL_RGB, GL_FLOAT, pixelcolorbuf);
  int U=(int)sqrt((float)themesh->faces.size()); 
  for (i = 0; i < width; i++)
	  for (j = 0; j < height; j++)  
		  if (pixelcolorbuf[i][j][2] == 0.0)   //pixel-level saliency 
		  {
			  index = int(U*(U-1)*pixelcolorbuf[i][j][0]+(U-1)*pixelcolorbuf[i][j][1]);
			  area[index]=area[index]+1;			
		  }

  double totalarea=0;
  for(i=0;i<themesh->faces.size();i++)
  {
  	if(area[i]!=0)
		m++;
	totalarea=totalarea+area[i];
  }
  
  //total area
  //S=width*height;  
  //Use the bounding circle
  //S=PI*(width/4.0)*(width/4.0); //  
  //S=(width/2.0)*(height/2.0)
  //S=GW*GH/4;
  //S=PI*25*25;
  
  S=(double) (GW*GH);//S=totalarea;  
    //S=PI*pow(GW/2.0,2.0);
  double area0=S-totalarea;
  entropy=-(area0/S)*(log(area0/S)/log(2.0));  
   for(i=0;i<themesh->faces.size();i++)
   	if (area[i]!=0)
		entropy=entropy-(area[i]/S)*(log(area[i]/S)/log(2.0));

  //normalization
  entropy=entropy/log(m+1)*log(2.0);
  delete [] area;
  return entropy;
}

