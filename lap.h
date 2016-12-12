#ifndef LAP_H
#define LAP_H
/************************************************************************
*
*  lap.h
   version 1.0 - 21 june 1996
   author  Roy Jonker, MagicLogic Optimization Inc.
   
   header file for LAP
*
**************************************************************************/

/*************** CONSTANTS  *******************/

#define BIG 10000000000000
#define m1 1000
#define n1 1000

/*************** TYPES      *******************/

  typedef int row;
  typedef int col;
  typedef int cost;

/*************** FUNCTIONS  *******************/

/*extern cost lap(int dim, //double [m][n],//double **assigncost,
               int *rowsol, int *colsol, double *u, double *v);*/
extern cost lblap(int dim, 
        //double assigncost[m][n],// cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v);

extern float lblapMB(int dim, 
        //double assigncost[m][n],// cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v);

extern void checklap(int dim,// double [m][n],//cost **assigncost,
                     int *rowsol, int *colsol, double *u, double *v);
//void ComputeShapeContext(double *F,double M1,double M2);
//void ComputeShapeContext(cost **F);
//void ComputeShapeContext(double assigncost[m][n]);
void ComputeShapeContext(int scale);
void ShapeContext() ; //Compute the original shape context matrix
void HistDistri(int scale);
#endif