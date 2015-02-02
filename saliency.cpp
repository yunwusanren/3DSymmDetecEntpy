#include "stdafx.h"
#include "mesh.hpp"
#include <vector>
#include <cmath>
#include <memory>
#include <d3dx10.h>
#include <algorithm>
#include <boost/array.hpp>
#include <functional>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <windows.h>
#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/foreach.hpp>
#include <boost/range.hpp>
#include "timestamp.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include "XForm.h"

using namespace std;
using namespace boost;
namespace ublas = boost::numeric::ublas;
#define foreach BOOST_FOREACH
int MainWindow;
Mesh themesh;
struct octtree_node;

struct octtree {
  vector<int> vert;
  auto_ptr<octtree_node> root;

  template<typename Functor>
  void search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance); // f(vert_type& v, float distance)
};

struct octtree_node {
  vector<int>::iterator from, to;
  D3DXVECTOR3 v0, v1;
  D3DXVECTOR3 center;
  float diagonal;
  array<auto_ptr<octtree_node>, 8> children; // arary of auto_ptr. IS THIS OK?

  template<typename Functor>
  void search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance);
};

template<float _D3DVECTOR::* memptr>
struct veccoord_cmp : unary_function<bool, int> {
  veccoord_cmp(Mesh& mesh, float threshold) : m(mesh), t(threshold) {}
  Mesh& m;
  float t;
  bool operator()(int i) const {
    const D3DXVECTOR3& v = m.vert_[i].p;
    return v.*memptr < t;
  }
};

template<typename Iterator>
auto_ptr<octtree_node> create_octtree_node(
  Mesh& m, Iterator from, Iterator to,
  const D3DXVECTOR3& v0, const D3DXVECTOR3& v1, int depth)
{
  // typeof(*from) == typeof(int), index of Mesh::vert_.
  auto_ptr<octtree_node> n(new octtree_node);
  n->from = from;
  n->to   = to;
  n->v0   = v0;
  n->v1   = v1;
  n->center   = (v0 + v1) / 2;
  D3DXVECTOR3 diag(v1 - v0);
  n->diagonal = D3DXVec3Length(&diag);
  if(/*depth != 0 && */to - from >= 16) {
    float vz[3];
    Iterator zit[3];
    vz[0] = v0.z; vz[2] = v1.z;
    vz[1] = (v0.z + v1.z) / 2;
    zit[0] = from; zit[2] = to;
    zit[1] = partition(from, to, veccoord_cmp<&D3DXVECTOR3::z>(m, vz[1]));
    for(int z = 0; z < 2; ++z) {
      float vy[3];
      Iterator yit[3];
      vy[0] = v0.y; vy[2] = v1.y;
      vy[1] = (v0.y + v1.y) / 2;
      yit[0] = zit[z]; yit[2] = zit[z+1];
      yit[1] = partition(zit[z], zit[z+1], veccoord_cmp<&D3DXVECTOR3::y>(m, vy[1]));
      for(int y = 0; y < 2; ++y) {
        float vx[3];
        Iterator xit[3];
        vx[0] = v0.x; vx[2] = v1.x;
        vx[1] = (v0.x + v1.x) / 2;
        xit[0] = yit[y]; xit[2] = yit[y+1];
        xit[1] = partition(yit[y], yit[y+1], veccoord_cmp<&D3DXVECTOR3::x>(m, vx[1]));
        for(int x = 0; x < 2; ++x) {
          D3DXVECTOR3 newv0(vx[x  ], vy[y  ], vz[z  ]);
          D3DXVECTOR3 newv1(vx[x+1], vy[y+1], vz[z+1]);
          n->children[z*4+y*2+x] = create_octtree_node(m, xit[x], xit[x+1], newv0, newv1, depth-1);
        }
      }
    }
  }
  return n;
}

auto_ptr<octtree> create_octtree(Mesh& m, const D3DXVECTOR3& v0, const D3DXVECTOR3& v1) {
  auto_ptr<octtree> o(new octtree);
  o->vert.resize(m.vert_.size());
  for(int i = 0; i < o->vert.size(); ++i) {
    o->vert[i] = i;
  }
  int depth = (int)(log((double)m.vert_.size()) - log(8.0));
  depth += 0;
  o->root = create_octtree_node(m, o->vert.begin(), o->vert.end(), v0, v1, depth);
  return o;
}

template<typename Functor>
void octtree_node::search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance) {
  D3DXVECTOR3 d3(center - this->center);
  const float d = D3DXVec3Length(&d3);
  if(d + diagonal < distance) {
    // all included
    foreach(int i, make_iterator_range(from, to)) {
      Mesh::vert_type& v = m.vert_[i];
      D3DXVECTOR3 vtoc(center - v.p);
      f(v, D3DXVec3Length(&vtoc));
    }
  }
  else if(d - diagonal > distance) {
    // not included
  }
  else {
    // may included
    if(children[0].get()) {
      // if I have children, give over.
      for(int i = 0; i < 8; ++i) {
        children[i]->search(m, f, center, distance);
      }
    }
    else {
      // if not, calc on each verts.
      foreach(int i, make_iterator_range(from, to)) {
        Mesh::vert_type& v = m.vert_[i];
        D3DXVECTOR3 vtoc(center - v.p);
        float d = D3DXVec3Length(&vtoc);
        if(d < distance) { f(v, d); }
      }
    }
  }
}

template<typename Functor>
void octtree::search(Mesh& m, Functor f, const D3DXVECTOR3& center, float distance) {
  root->search(m, f, center, distance);
}

struct sal_func {
  sal_func(float epsilon, double* sigma, double* sigma_denom)
       : eps(epsilon), sig(sigma), sigd(sigma_denom)
    {
      for(int i = 0; i < 13; ++i) {
        exp_rcp2_s[i] = exp(1.f / (i*i));
      }
    }
  double *sig, *sigd;
  float eps;
  float exp_rcp2_s[13];
  void operator()(Mesh::vert_type& vj, float dist) const {
    if(dist == 0.f) { return; }
    float distance2 = dist*dist;
    int distancei = (int)(fabs(dist / eps));
    float e = exp(-distance2 / (eps*eps) / 2);

    for(int s = distancei; s < 13; ++s) {
      float e2 = e * exp_rcp2_s[s];
      float curv_mean = (vj.curv_tensor[0] + vj.curv_tensor[2]) / 2;
      sig[s]  += e2 * curv_mean;
      sigd[s] += e2;
    }
  }
};

void Mesh::calc_saliency() {
  // make bounding-box
  float epsilon;
  D3DXVECTOR3 bound_min( FLT_MAX,  FLT_MAX,  FLT_MAX);
  D3DXVECTOR3 bound_max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
  {
    for(int i = 0; i < vert_.size(); ++i) {
      vert_type& v = vert_[i];
      D3DXVec3Minimize(&bound_min, &bound_min, &v.p);
      D3DXVec3Maximize(&bound_max, &bound_max, &v.p);
    }
    D3DXVECTOR3 diagonal(bound_max-bound_min);
    epsilon = D3DXVec3Length(&diagonal) * 0.003;
  }
  auto_ptr<octtree> ot(create_octtree(*this, bound_min, bound_max));

  array<float, 5> zeros = { };
  vector<array<float, 5> > saliencies(vert_.size(), zeros);

  // compute saliencies on each scale level.
//#pragma omp parallel for
  for(int i = 0; i < vert_.size(); ++i) {
    // print progress
    if(i % (vert_.size() / 10) == 0) {
      cout << "saliency " << (i / (vert_.size() / 10)) << "0%" << endl;
    }

    vert_type& vi = vert_[i];
    array<double, 13> sigma = {};
    array<double, 13> sigma_denom = {};
    sal_func f(epsilon, sigma.c_array(), sigma_denom.c_array());
    ot->search(*this, f, vi.p, epsilon * 12);

    for(int s = 0; s < 13; ++s) {
      if(sigma_denom[s] == 0.0) { sigma_denom[s] += 1.0; }
      sigma[s] /= sigma_denom[s];
    }

    for(int e = 2; e <= 6; ++e) {
      saliencies[i][e-2] = fabs(sigma[e] - sigma[e*2]);
    }
  }

  // normalize saliencies
  vector<float> max_saliencies(5, 0.0);
  for(int i = 0; i < vert_.size(); ++i) {
    for(int e = 0; e < 5; ++e) {
      max_saliencies[e] = max(max_saliencies[e], saliencies[i][e]);
    }
  }
  for(int i = 0; i < vert_.size(); ++i) {
    for(int e = 0; e < 5; ++e) {
      saliencies[i][e] /= max_saliencies[e];
    }
  }

  // find local maximams
  array<float, 5> local_max = {};
  array<float, 5> local_max_denom = {};
  for(int i = 0; i < vert_.size(); ++i) {
    float maxs[] = { -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX, -FLT_MAX };
    for(int j = 0; j < adjvert_[i].size(); ++j) {
      for(int e = 0; e < 5; ++e) {
        maxs[e] = max(maxs[e], saliencies[j][e]);
      }
    }

    for(int e = 0; e < 5; ++e) {
      float sal = saliencies[i][e];
      if(maxs[e] == sal && sal != 1.0) {
        local_max[e] += sal;
        local_max_denom[e] += 1.f;
      }
    }
  }
  for(int e = 0; e < 5; ++e) {
    if(local_max_denom[e] == 0.0) { local_max_denom[e] = 1.0; }
    local_max[e] /= local_max_denom[e];
  }

  // compute final saliency
  float scale_coeff[5];
  for(int e = 0; e < 5; ++e) {
    scale_coeff[e] = pow(1.0 - local_max[e], 2);
  }
  vector<float> total_saliencies(vert_.size(), 0.f);
  for(int i = 0; i < vert_.size(); ++i) {
    for(int e = 0; e < 5; ++e) {
      total_saliencies[i] += saliencies[i][e] * scale_coeff[e];
    }
  }

  foreach(float& sal, total_saliencies) {
    if(!(sal <= 0.f || 0.f <= sal)) {
      // v.saliency is NaN
      sal = 1.f;
    }
  }

  // smooth
  for(int i = 0; i < vert_.size(); ++i) {
    vert_type& v = vert_[i];
    v.saliency = total_saliencies[i];
    for(int j = 0; j < adjvert_[i].size(); ++j) {
      v.saliency += total_saliencies[adjvert_[i][j]];
    }
    v.saliency /= adjvert_[i].size() + 1;
    printf("%f ",v.saliency);
  }
}

#define DIM 20   
#define TWOPI 6.2831852
xform xf;
//GLCamera camera;
 float ph=(1+sqrt(5.0))/2;  
float ViewPort[20][3]={{1,1,1},{1,1,-1},{1,-1,1},{1,-1, -1},{-1,1,1},{-1, 1,-1},{-1, -1,1},{-1,-1,-1},
 {0,1/ph,ph},{0,1/ph, -ph},{0,-1/ph,ph},{0, -1/ph, -ph},{1/ph,ph,0},{1/ph,-ph,0},{-1/ph,ph,0},
{-1/ph,-ph,0},{ph,0,1/ph},{ph,0,-1/ph},{-ph, 0,1/ph},{-ph,0,-1/ph}};     
  float axis1[3]={1,0,0};  //rotate about X axis for phai
  float axis2[3]={0,1,0};  //rotate about Y axis for theta
  float axis3[3]={0,0,1};  //rotate about Z axis for theta
void dispmodel()
{	 
	int i,num1,num2; 

	//Setup the camera
	glMatrixMode(GL_MODELVIEW); 
	glLoadIdentity(); 
	gluLookAt(0,0,20,0,0,0,0,1,0); 
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);	//glShadeModel(GL_FLAT);	
	glShadeModel(GL_SMOOTH);
	glPolygonMode(GL_FRONT, GL_FILL);
	glPolygonMode(GL_BACK, GL_FILL);

	glPointSize(1);
	/*
	for (i=0;i<themesh.vert_.size();i++)
	{
		num1=i/256;
		num2=i-256*num1;
		glColor3f(num2/255.0,num1/255.0,0.0);
	
		glBegin(GL_POINTS);
		glVertex3f(themesh.vert_[i].p.x,themesh.vert_[i].p.y,themesh.vert_[i].p.z);
		glEnd();
	} 
	*/
	for (i=0;i<(int)themesh.face_.size();i++)
	{
		glPointSize(1);
		glBegin(GL_TRIANGLE_STRIP);

		num1=themesh.face_[i].vpos[0]/256;
		num2=themesh.face_[i].vpos[0]-256*num1;
		glColor3f(num2/255.0,num1/255.0,0.0);
		D3DXVECTOR3 v;
		v=xf*themesh.vert_[themesh.face_[i].vpos[0]].p;
		glNormal3f(themesh.vert_[themesh.face_[i].vpos[0]].n.x,themesh.vert_[themesh.face_[i].vpos[0]].n.y,themesh.vert_[themesh.face_[i].vpos[0]].n.z);
		glVertex3f(v.x,v.y,v.z);	

		num1=themesh.face_[i].vpos[1]/256;
		num2=themesh.face_[i].vpos[1]-256*num1;
		glColor3f(num2/255.0,num1/255.0,0.0);
		v=xf*themesh.vert_[themesh.face_[i].vpos[1]].p;
		glNormal3f(themesh.vert_[themesh.face_[i].vpos[1]].n.x,themesh.vert_[themesh.face_[i].vpos[1]].n.y,themesh.vert_[themesh.face_[i].vpos[1]].n.z);
		glVertex3f(v.x,v.y,v.z);	

		num1=themesh.face_[i].vpos[2]/256;
		num2=themesh.face_[i].vpos[2]-256*num1;
		glColor3f(num2/255.0,num1/255.0,0.0);
		v=xf*themesh.vert_[themesh.face_[i].vpos[2]].p;
		glNormal3f(themesh.vert_[themesh.face_[i].vpos[2]].n.x,themesh.vert_[themesh.face_[i].vpos[2]].n.y,themesh.vert_[themesh.face_[i].vpos[2]].n.z);
		glVertex3f(v.x,v.y,v.z);	
		glEnd();
	} 
}
/*---------------------------------------------------------------------------
  Resize the window
  ---------------------------------------------------------------------------*/
void reshape(int w, int h)
{
  //GW = w;
 // GH = h;
  //AR = (float)(w)/(float)(h);
  glViewport(0, 0, w, h);               /* Set Viewport */   
  glutPostRedisplay();
}
void disp()
{
   glEnable(GL_DEPTH_TEST); 
   // Just clean the screen
   glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT); 
   
   // lighting
   GLfloat light_position[]={ 1.0,1.0,1.0,0.0 };

   //GLfloat light_position[]={ 0.0,0.0,10.0,0.0 };
   GLfloat white_light[]={ 1.0,1.0,1.0,1.0 };
   GLfloat lmodel_ambient[]={ 1.0,0.1,0.1,1.0 };

   glLightfv(GL_LIGHT0, GL_POSITION, light_position);
   glLightfv(GL_LIGHT0, GL_DIFFUSE, white_light);
   glLightfv(GL_LIGHT0, GL_SPECULAR, white_light);

   glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
   glEnable(GL_LIGHTING);glEnable(GL_LIGHT0);
   
   glMatrixMode(GL_PROJECTION); 
   glLoadIdentity(); 
   gluLookAt(0,0,20,0,0,0,0,1,0); 

   glMatrixMode(GL_PROJECTION); 
   glLoadIdentity(); 
   glOrtho(-2, 2, -2, 2, -100, 100);
   
   glMatrixMode(GL_MODELVIEW); 
   glLoadIdentity(); 
   glClearColor(1.0,1.0, 1.0,1.0);
   glDisable(GL_LIGHTING);
   dispmodel(); 
   glutSwapBuffers();
   glutPostRedisplay();
}
 
void resetview()
{
	//xf=xform::trans(-themesh->bsphere.center);
	//xf = xf*xform::trans(0, 0, -10.0f * themesh->bsphere.r);
	//camera.stopspin();
	xform PCA(1,0,0,0,0,1,0,0,0,0,1,0,0, 0, 0, 1);
	xf=PCA;
}

void ComputeCurrentCubicView(int i)
{
	float phai1,theta1,x,y,z;
	
	x=ViewPort[i][0];y=ViewPort[i][1];z=ViewPort[i][2];
	phai1=-acos(y/sqrt(x*x+y*y+z*z));
	if (phai1 < 0)
		phai1 = TWOPI - fabs(phai1);
	theta1=-atan2(z, x);
	if (theta1 < 0)
		theta1 = TWOPI - fabs(theta1);
	resetview();			
	xf=xf*xform::rot(theta1, axis2[0], axis2[1], axis2[2])*xform::rot(phai1, axis1[0], axis1[1], axis1[2]);
	glutSetWindow( MainWindow);
	disp(); //glutPostRedisplay(); 	
}

float pixelcolorbuf[500][500][3] ={0}; 
float calVisibleSalency()
{
  int i, j, k,index;
  float visiblesaliency=0.0;
  int View[4];
  glutSetWindow( MainWindow);
  glGetIntegerv(GL_VIEWPORT, View);
  int width = View[2], height = View[3];

  for (i = 0; i < width; i++)
    for (j = 0; j < height; j++)
      for (k = 0; k < 3; k++)
		  pixelcolorbuf[i][j][k] = 0;
  
  glReadBuffer(GL_FRONT);
  glReadPixels(View[0], View[1], width, height, GL_RGB, GL_FLOAT, pixelcolorbuf);
  for (i = 0; i < width; i++)
	  for (j = 0; j < height; j++)  
		  //R or G component are not zero
		 // if ((pixelcolorbuf[i][j][0] != 0 || pixelcolorbuf[i][j][1] != 0) && (pixelcolorbuf[i][j][2] == 0.0))
		  if (pixelcolorbuf[i][j][2] == 0.0)   //pixel-level saliency 
		  {
			  index = int(255 *pixelcolorbuf[i][j][1]) *256+int(255*pixelcolorbuf[i][j][0]);
			  if (index>= themesh.vert_.size())
			  	;//printf("out of range!\n");
			  else
			  	visiblesaliency=visiblesaliency+themesh.vert_[index].saliency;
		  }
  return visiblesaliency;
}
float maxVisibleSaliency()
{
	int i,pos;
	float maxS=0,s;
	for(i=0;i<DIM;i++)
	{
		ComputeCurrentCubicView(i);
		s=calVisibleSalency();
		printf("Pose:%d,Saliency:%f\n",i,s);
		if (maxS<s)
		{
		   maxS=s;
		   pos=i;
		 }
	}	
	ComputeCurrentCubicView(pos);
	printf("Optimal pose: %f, %f, %f\n",ViewPort[pos][0],ViewPort[pos][1],ViewPort[pos][2]);
	return maxS;
}
//Load A Model file
int  modelfile()
{

	OPENFILENAME OpenFileName;
	char szFile[MAX_PATH];

    //Create File Browser// 
	szFile[0] = 0;
	
    HWND hwnd= GetActiveWindow();              // owner window
    //printf("\nhwnd=%d\n",hwnd);
	OpenFileName.lStructSize = sizeof( OPENFILENAME );
	OpenFileName.hwndOwner = (HWND)hwnd;
	//OpenFileName.lpstrFilter =(LPCSTR)"\PLY(*.ply)\0*.ply\0All Files(.*)\0*.*\0";
	//OpenFileName.lpstrFilter =(LPCSTR)"\OFF(*.off)\0*.off\0All Files(.*)\0*.*\0";
	OpenFileName.lpstrFilter =(LPCSTR)"\OBJ(*.obj)\0*.obj\0All Files(.*)\0*.*\0";
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
		printf("\n%s\n",szFile);	
		//delete themesh;
		filebuf fb;
		fb.open (szFile,ios::in); //armadillo.obj
       	istream ifs(&fb);
		timestamp t1 = now();
		//Compute mesh saliency 
		themesh.load_obj(ifs);
		themesh.make_adjinfo();
		themesh.calc_normal();
		themesh.make_ntb();
		themesh.calc_curv_dir();
		themesh.smooth_curv();
		themesh.calc_saliency();	
		printf("\rMesh Saliency Computation Elapsed time: %.2f sec.\n",  (now() - t1));
		return FALSE;
	}
	else 
	{   
		printf("Loading file  Failed!");
		return FALSE;
	}	
}

/*---------------------------------------------------------------------------
Mouse Press p callback routine
--------------------------------------------------------------------------*/
void mymouse(int button, int state, int x, int y)
{
  int i;
  timestamp t1 = now();
  if (glutGetWindow() == MainWindow)
  {
	  if (state == GLUT_DOWN && button == GLUT_MIDDLE_BUTTON)
	  {

		  float totalsaliency=0.0,visiblesaliency;
/*		  //Compute the total mesh saliency of all points
		  for(i=0;i<themesh.vert_.size();i++)
			  totalsaliency=totalsaliency+themesh.vert_[i].saliency;
		  printf("Total Saliency=%f\n",totalsaliency);
*/		  
		  //Compute the total mesh saliency of visible points
		  //visiblesaliency=calVisibleSalency();
		  visiblesaliency=maxVisibleSaliency();
		  printf("Visible Saliency=%f\n",visiblesaliency);
		  
		  printf("\rMaximal Visible Saliency Computation Elapsed time: %.2f sec.\n",  (now() - t1));
	  }
	  else if (state == GLUT_DOWN && button == GLUT_RIGHT_BUTTON)
	  {
	     modelfile();
	  }
  }   
}

/*---------------------------------------------------------------------------
Resize the window
---------------------------------------------------------------------------*/
void reshape1(int w, int h)
{
  //GW = w;
  //GH = h;
  //AR = (float)(w) / (float)(h);
  glViewport(0, 0, w, h); /* Set Viewport */
  glutPostRedisplay();
}

int main()
{
	filebuf fb;


	//fb.open ("D:\\Programming\\Saliency\\saliency\\Debug\\bunny.obj",ios::in);
	//fb.open ("D:\\Programming\\Saliency\\saliency\\Debug\\backdoor.obj",ios::in);
	//fb.open ("D:\\Programming\\Saliency\\saliency\\Debug\\horse.obj",ios::in);
	fb.open ("D:\\Programming\\Saliency\\saliency\\Debug\\santa.obj",ios::in); //armadillo.obj

	istream ifs(&fb);

	timestamp t1 = now();
	//Compute mesh saliency 
	themesh.load_obj(ifs);
	themesh.make_adjinfo();
	themesh.calc_normal();
	themesh.make_ntb();
	themesh.calc_curv_dir();
	themesh.smooth_curv();
	themesh.calc_saliency();
	printf("\rMesh Saliency Computation Elapsed time: %.2f sec.\n",  (now() - t1));
	
	//initialize the GLUT
	glutInitWindowSize(500,500);            	//glutInitWindowSize(760, 700);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	//glutInit(&argc, argv);
	char windowname[30]="Salency-based alignment";
	glutInitWindowPosition(50,100);
	MainWindow=glutCreateWindow(windowname);
	glViewport(0,0,500,500 );
	glClearColor(1.0,1.0,1.0,1.0);

	glutDisplayFunc(disp);
	glutReshapeFunc(reshape);
	glutMouseFunc(mymouse);
	//glutMotionFunc(mymotion); 
	//glutKeyboardFunc(KeyFunc);
	glutMainLoop();	
	
}
