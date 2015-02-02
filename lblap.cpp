//************************************************************************
//
//  lap.cpp
//   version 1.0 - 4 September 1996
//   author: Roy Jonker @ MagicLogic Optimization Inc.
//   e-mail: roy_jonker@magiclogic.com

//   Code for Linear Assignment Problem, according to 
   
//   "A Shortest Augmenting Path Algorithm for Dense and Sparse Linear   
//    Assignment Problems," Computing 38, 325-340, 1987
   
//   by
   
//   R. Jonker and A. Volgenant, University of Amsterdam.
//
//*************************************************************************/

//global assigncost[m][n];

//function lapcost=lap(dim,assigncost,rowsol,colsol,u,v)
//float lblap(int dim,double assigncost[m][n])
#include "stdafx.h"
#include "system.h"
#include "gnrl.h"
#include "lap.h"
#include "GL/glui.h"
#include "Wire.h"
#include <climits>

float  flapcost=0;
extern GLUI *glui;
extern cost assigncost[2563][2563]; //[1282][1282];
extern int startdim,enddim;
#define NUM 21  // 200  500      //Sampling points

void HistDistri(int scale);
float lblapMB(int dim, 
        //double assigncost[m][n],// cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v)
{
  int unassignedfound;
  row i, imin, numfree = 0, prvnumfree, f, i0, k, freerow;
  col j, j1, j2, endofpath, last, low, up;
  cost mini, h, umin, usubmin, v2, lapcost;

  //row rowsol[NUM + 1];
  //col colsol[NUM + 1];
  //cost u[NUM + 1];
  //cost v[NUM + 1];

  row unused[NUM + 1];   // list of unassigned rows
  col collist[NUM + 1];  // list of columns to be scanned in various ways
  cost d[NUM + 1];       // 'cost-distance' in augmenting path calculation
  row pred[NUM + 1];     // row-predecessor of column in augmenting/alternating path

  for (i = 1; i <= NUM; i++) {
    rowsol[i] = 0;
  }
  // COLUMN REDUCTION 
  for (j = NUM; j >= 1; j--) {  // reverse order gives better results
    // find minimum cost over rows
    mini = assigncost[1][j];
    imin = 1;
    for (i = 2; i <= NUM; i++) {
      if (assigncost[i][j] < mini) {
        mini = assigncost[i][j];
        imin = i;
      }
    }
    v[j] = mini;

    if (rowsol[imin] == 0) {  // init assignment if minimum row assigned for first time
      rowsol[imin] = j;
      colsol[j] = imin;
    } else {
      if (rowsol[i] > 0) {
        rowsol[i] = -rowsol[i];
      }
      colsol[j] = 0;  // row already assigned, column not assigned
    }
  }

  numfree = 0;
  // REDUCTION TRANSFER
  for (i = 1; i <= NUM; i++) {
    if (rowsol[i] == 0) {  // fill list of unassigned 'free' rows
      numfree++;
      unused[numfree] = i;
    }
    if (rowsol[i] < 0)
      rowsol[i] = -rowsol[i];
    else {  // transfer reduction from rows that are assigned once
      j1 = rowsol[i];
      //mini = BIG;  // error
      mini = INT_MAX;  // MB fix
      for (j = 1; j <= NUM; j++) {
        if (j != j1)
          if ((assigncost[i][j] - v[j]) < mini)
            mini = assigncost[i][j] - v[j];
      }
      v[j1] -= mini;
    }
  }

  // AUGMENTING ROW REDUCTION 
  int loopcnt = 0;  // loop to be done twice
  while (loopcnt <= 1) {
    // scan all free rows
    // in some cases, a free row may be replaced with another one to be scanned next
    k = 1;
    prvnumfree = numfree;
    numfree = 0;  // start list of rows still free after augmenting row reduction
    while (k <= prvnumfree) {
      i = unused[k];
      k++;
      umin = assigncost[i][1] - v[1];
      j1 = 1;
      //usubmin = BIG;  // error
      usubmin = INT_MAX;  // MB fix
      // find minimum and second minimum reduced cost over columns
      for (j = 2; j <= NUM; j++) {
        h = assigncost[i][j] - v[j];
        if (h < usubmin) {
          if (h >= umin) {
            usubmin = h;
            j2 = j;
          } else {
            usubmin = umin;
            umin = h;
            j2 = j1;
            j1 = j;
          }
        }
      }

      i0 = colsol[j1];

      if (umin < usubmin) {
        // change the reduction of the minimum column to increase the minimum
        // reduced cost in the row to the subminimum
        v[j1] += umin - usubmin;
      } else if (i0 > 0) {  // minimum column j1 is assigned
        // swap columns j1 and j2, as j2 may be unassigned
        j1 = j2;
        i0 = colsol[j1];
      }

      if (i0 > 0) {  // minimum column j1 assigned earlier
        if (umin < usubmin) {
          // put in current k, and go back to that k
          // continue augmenting path i - j1 with i1
          k--;
          unused[k] = i0;
        } else {
          // no further augmenting reduction possible
          // store i1 in list of free rows for next phase
          numfree++;
          unused[numfree] = i0;
        }
      }
      rowsol[i] = j1;
      colsol[j1] = i;
    }   //while (k < prvnumfree)
    loopcnt++;
  }   //while (loopcnt <=2)   

  for (j = 1; j <= NUM; j++) {
    collist[j] = j;
  }
  // AUGMENT SOLUTION for each free row
  prvnumfree = numfree;
  for (f = 1; f <= prvnumfree; f++) {
    freerow = unused[f];  // start row of augmenting path
    low = 1;  // columns in 0..low-1 are ready, now none
    up = 1;  // columns in low..up-1 are to be scanned for current minimum, now none

    // Dijkstra shortest path algorithm
    // runs until unassigned column added to shortest path tree
    for (j = 1; j <= NUM; j++) {
      d[j] = assigncost[freerow][j] - v[j];
      pred[j] = freerow;
    }

    // columns in up..NUM-1 are to be considered later to find new minimum, 
    // at this stage the list simply contains all columns
    unassignedfound = 0;

    while (!unassignedfound) {
      if (up == low) {  // no more columns to be scanned for current minimum
        last = low - 1;
        // scan columns for up..NUM-1 to find all indices for which new minimum occurs
        // store these indices between low..up-1 (increasing up)              
        mini = d[collist[up]];
        up++;
        for (k = up; k <= NUM; k++) {
          j = collist[k];
          h = d[j];
          if (h <= mini) {
            if (h < mini) {  // new minimum 
              up = low;  // restart list at index low
              mini = h;
            }
            // new index with same minimum, put on undex up, and extend list
            collist[k] = collist[up];
            collist[up] = j;
            up++;
          }
        }

        // check if any of the minimum columns happens to be unassigned
        // if so, we have an augmenting path right away
        for (k = low; k < up; k++) {
          endofpath = collist[k];
          if (colsol[collist[k]] == 0) {
            unassignedfound = 1;
            break;
          }
        }
      }   //if (up==low)

      if (!unassignedfound) {
        // update 'distances' between freerow and all unscanned columns, via next scanned column
        j1 = collist[low];
        low++;
        i = colsol[j1];
        h = assigncost[i][j1] - v[j1] - mini;

        for (k = up; k <= NUM; k++) {
          j = collist[k];
          endofpath = j;
          v2 = assigncost[i][j] - v[j] - h;
          if (v2 < d[j]) {
            d[j] = v2;
            pred[j] = i;
            if (v2 == mini) {  // new column found at same minimum value
              if (colsol[j] == 0) {
                // if unassigned, shortest augmenting path is complete                             
                unassignedfound = 1;
                break;
                // else add to list to be scanned right away
              } else {
                collist[k] = collist[up];
                collist[up] = j;
                up++;
              }
            }
          }
        }   //(for k = up:NUM)
      }   //if (~unassignedfound)
    }   //while (~unassignedfound)

    // update column prices
    for (k = 1; k <= last; k++) {
      j1 = collist[k];
      v[j1] += d[j1] - mini;
    }

    // reset row and column assignments along the alternating path
    i = pred[endofpath];
    colsol[endofpath] = i;
    j1 = endofpath;
    endofpath = rowsol[i];
    rowsol[i] = j1;

    while (i != freerow) {
      i = pred[endofpath];
      colsol[endofpath] = i;
      j1 = endofpath;
      endofpath = rowsol[i];
      rowsol[i] = j1;
    }
  }   //for f = 1:numfree 

  // calculate optimal cost
  lapcost = 0;
  int nmatches = 0;  //number of successful matches
  for (i = 1; i <= NUM; i++) {
    j = rowsol[i];
    if (j > 0) {
      u[i] = assigncost[i][j] - v[j];
      lapcost += assigncost[i][j];
      nmatches++;
    }
  }

  return float(lapcost) / Examplifier;  //to show the result
}
  
cost lblap(int dim, 
        //double assigncost[m][n],// cost **assigncost,
        col *rowsol, 
        row *colsol, 
        cost *u, 
        cost *v)
{
//comm1
//assigncost=F;
//dim1=max(size(F));
//dim=min(size(F));
//assigncost=F;
// input:
// dim        - problem size
// assigncost - cost matrix

// output:
// rowsol     - column assigned to row in solution
// colsol     - row assigned to column in solution
// u          - dual variables, row reduction numbers
// v          - dual variables, column reduction numbers

  //int unassignedfound;
  //row  i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
  //col  j, j1, j2, endofpath, last, low, up, *collist, *matches;
  //cost min, h, umin, usubmin, v2, *d;
  int unassignedfound;
  row  i, imin, numfree = 0, prvnumfree, f, i0, k, freerow, *pred, *free;
  col  j, j1, j2, endofpath, last, low, up, *collist;
  cost mini,h, umin, usubmin, v2, *d,lapcost;


  free = new row[dim+1];       // list of unassigned rows.
  collist = new col[dim+1];    // list of columns to be scanned in various ways.
  d = new cost[dim+1];         // 'cost-distance' in augmenting path calculation.
  pred = new row[dim+1];       // row-predecessor of column in augmenting/alternating path.

  for(i=1;i<=dim;i++)
  {
	  free[i] =0;       // list of unassigned rows.
      collist[i] =0;    // list of columns to be scanned in various ways.
	  d[i] =0;        // 'cost-distance' in augmenting path calculation.
	  pred[i] =0;      // row-predecessor of column in augmenting/alternating path.
	  rowsol[i] =0;
	  colsol[i] =0;
	  u[i] =0;
	  v[i] =0;
  }
  
  // COLUMN REDUCTION 
  for(j=dim;j>=1;j--) // reverse order gives better results.
  {
	  // find minimum cost over rows.
      collist[j]=j;
      mini = assigncost[1][j];     
      imin = 1;
      for(i=2;i<=dim;i++) 
		  if (assigncost[i][j] < mini) 
		  {
			  mini = assigncost[i][j];        
              imin = i;
		  }     
      v[j] = mini;       

      if (rowsol[imin] == 0)
	  { // init assignment if minimum row assigned for first time.
          rowsol[imin] = j; 
          colsol[j] = imin; 
	  }
      else
	  {
		  rowsol[i] = -abs(rowsol[i]);  
          colsol[j] =0;        // row already assigned, column not assigned.
	  } 
  }
 numfree=0;
  // REDUCTION TRANSFER
  for(i=1;i<=dim;i++)
  {
	  if (rowsol[i] ==0)     // fill list of unassigned 'free' rows.
	  {
		  numfree=numfree+1;
          free[numfree] = i;
	  }
      if (rowsol[i]<0)   
          rowsol[i]=-rowsol[i];
      else
	  {   // transfer reduction from rows that are assigned once.
          j1 = rowsol[i];mini = BIG;
           for(j=1;j<=dim;j++)
		   { 
			   if (j != j1)
		    if ((assigncost[i][j] - v[j])<mini) 
				mini = assigncost[i][j]- v[j];
		   }
          v[j]=v[j]-mini;
     }
  }
  
  // AUGMENTING ROW REDUCTION 
  int loopcnt = 0;           // loop to be done twice.
  while (loopcnt <=1)   
  {          
      // scan all free rows.
      // in some cases, a free row may be replaced with another one to be scanned next.
      k = 1; prvnumfree = numfree; numfree = 0; // start list of rows still free after augmenting row reduction.
      while (k <=prvnumfree )//&& umin ~= usubmin)
	  {
          i = free[k]; k=k+1;umin = assigncost[i][1] - v[1];j1 = 1;usubmin = BIG;       
          // find minimum and second minimum reduced cost over columns.
          for (j=2;j<=dim;j++) 
		  {
			  h = assigncost[i][j] - v[j];
              if (h<usubmin) 
			  {
                  if (h>=umin) 
				  {
					  usubmin = h; 
					  j2 = j;
				  }
                  else 
				  {
					  usubmin = umin; 
                      umin = h; 
                      j2 = j1; 
                      j1 = j;
				  }
			  }
		  } 
          
          i0 = colsol[j1];
          
          if (umin<usubmin) 
              // change the reduction of the minimum column to increase the minimum
              // reduced cost in the row to the subminimum.
              v[j1] = v[j1] + umin - usubmin;
          else
              if (i0 >0)         // minimum column j1 is assigned.
			  {
				  // swap columns j1 and j2, as j2 may be unassigned.
                  j1 = j2; 
                  i0 = colsol[j1];
			  }
                          
          if (i0 >0)           // minimum column j1 assigned earlier.
		  {
			  if (umin<usubmin) 
			  {
				  // put in current k, and go back to that k.
                  // continue augmenting path i - j1 with i1.
                  k=k-1;
				  free[k]= i0; 
			  }
              else 
			  {
				  // no further augmenting reduction possible.
                  // store i1 in list of free rows for next phase.
                  numfree=numfree+1;
                  free[numfree]= i0;     
			  }
		  }
          rowsol[i]=j1;colsol[j1]=i;
	  }  //while (k < prvnumfree)
      loopcnt=loopcnt+1;
 } //while (loopcnt <=2)   
  
  // AUGMENT SOLUTION for each free row.
  prvnumfree = numfree;  
  for (f =1;f<=prvnumfree;f++)
  {
      freerow = free[f];       // start row of augmenting path.
      low = 1; // columns in 0..low-1 are ready, now none.
      up = 1;  // columns in low..up-1 are to be scanned for current minimum, now none.
            
      // Dijkstra shortest path algorithm.
      // runs until unassigned column added to shortest path tree.
	  for (j = 1;j<=dim;j++)
	  {
          d[j] = assigncost[freerow][j] - v[j]; 
          pred[j] = freerow;
	  }     
            
      // columns in up..dim-1 are to be considered later to find new minimum, 
      // at this stage the list simply contains all columns
      //unassignedfound = FALSE;
      unassignedfound = 0;
      
      while (!unassignedfound )
	  {
		  if (up == low)         // no more columns to be scanned for current minimum.
		  {
			  last = low - 1;          
              // scan columns for up..dim-1 to find all indices for which new minimum occurs.
              // store these indices between low..up-1 (increasing up).              
              mini = d[collist[up]]; up=up+1;  //up,d(collist(up,1),1)
              for (k = up;k<=dim;k++)
			  {
				  j = collist[k]; 
                  h = d[j];
				  if (h<=mini)       
				  {
					  if (h<mini)   // new minimum. 
					  {
						  up = low;      // restart list at index low.
                          mini = h;
					  }                     
                      // new index with same minimum, put on undex up, and extend list.
                      collist[k] = collist[up]; 
                      collist[up] = j; 
                      up=up+1;
				  }
			  }
              
              // check if any of the minimum columns happens to be unassigned.
              // if so, we have an augmenting path right away.
              for (k = low;k<=up-1;k++)
			  {
				  endofpath = collist[k];
                  if (colsol[collist[k]]== 0)                      
				  {
					  unassignedfound = 1;
                      break;
				  }
			  }
		  }  //if (up==low)
          
          if (!unassignedfound) 
		  {
			  // update 'distances' between freerow and all unscanned columns, via next scanned column.
              j1 = collist[low]; 
              low=low+1; 
              i = colsol[j1]; 
              h = assigncost[i][j1] - v[j1] - mini; 
              
              for (k = up;k<=dim;k++)
			  {
				  j = collist[k]; 
                  endofpath = j;
                  v2 = assigncost[i][j] - v[j] - h;
                  if (v2 <d[j]) 
				  {
					  d[j] = v2;
                      pred[j] = i;
                      if (v2 == mini)   // new column found at same minimum value
					  {
						  if (colsol[j]==0) 
						  {
							  // if unassigned, shortest augmenting path is complete.                             
                              unassignedfound = 1;
                              break; 
                              // else add to list to be scanned right away.
						  }
                          else    
						  {
							  collist[k] = collist[up]; 
                              collist[up]= j; 
                              up=up+1;
						  }
					  }             
				  }
			  } //(for k = up:dim)
		  }//if (~unassignedfound)
	  } //while (~unassignedfound)
      
      // update column prices.
	  for (k =1;k<=last;k++)
	  {
          j1 = collist[k]; 
          v[j1] = v[j1] + d[j1] - mini;
	  }
      
      // reset row and column assignments along the alternating path.
      i = pred[endofpath]; 
      colsol[endofpath] = i; 
      j1 = endofpath; 
      endofpath = rowsol[i];
      rowsol[i] = j1;
      
	  while (i!= freerow)
	  {
          i = pred[endofpath]; 
          colsol[endofpath] = i; 
          j1 = endofpath; 
          endofpath = rowsol[i];
          rowsol[i] = j1;
	  }
  }//for f = 1:numfree 
  
  // calculate optimal cost.
  lapcost = 0;
  int nmatches=0;    //number of successful matches
  for (i =1;i<=dim;i++)   //for (i =0;i<startdim;i++) 
  {
	  j = rowsol[i];
	  if (j>0)
	  {u[i] = assigncost[i][j] - v[j];
	  lapcost = lapcost + assigncost[i][j];}
  }
  //lapcost=(float)lapcost/Examplifier;
  // free reserved memory.
  delete[] pred;
  delete[] free;
  delete[] collist;  
  delete[] d;
  //lapcost=int(lapcost/1000000/double(nmatches)*100);
  flapcost=float(lapcost)/Examplifier;  //to show the result
  lapcost=cost(flapcost);
  return lapcost;
}
  

void checklap(int dim, //cost **assigncost,
              col *rowsol, row *colsol, cost *u, cost *v)
{
  row  i;
  col  j;
  cost lapcost = 0, redcost = 0;
  int *matched;
//  char wait;
  
  matched = new int[dim];
  
  for (i = 1; i <=dim; i++)  
    for (j = 1; j <=dim; j++)  
      if ((redcost = assigncost[i][j] - u[i] - v[j]) < 0)
      {
        printf("\n");
        printf("negative reduced cost i= %d,j= %d,redcost %f\n", i, j, redcost);
        printf("\n\ndim %5d - press key\n", dim);
       // scanf("%d", &wait);
        break; 
      }

  for (i = 1; i <=dim; i++)  
    if ((redcost = assigncost[i][rowsol[i]] - u[i] - v[rowsol[i]]) != 0)
    {
      printf("\n");
      printf("non-null reduced cost i= %d,soli= %d,redcost= %f\n", i, rowsol[i], redcost);
      printf("\n\ndim %5d - press key\n", dim);
      //scanf("%d", &wait);
      break; 
    }
  
  for (j = 1; j <=dim; j++)  
    matched[j] = FALSE;
    
  for (i = 1; i <=dim; i++)  
    if (matched[rowsol[i]])
    {
      printf("\n");
      printf("column matched more than once - i= %d,soli= %d\n", i, rowsol[i]);
      printf("\n\ndim %5d - press key\n", dim);
     // scanf("%d", &wait);
      break; 
    }
    else
      matched[rowsol[i]] = TRUE;
      
    
  for (i = 1; i <=dim; i++)  
    if (colsol[rowsol[i]] != i)
    {
      printf("\n");
      printf("error in row solution i=%d,rowsol[i]= %d,colsol[rowsol[i]]= %d\n", i, rowsol[i], colsol[rowsol[i]]);
      printf("\n\ndim %5d - press key\n", dim);
      //scanf("%d", &wait);
      break; 
    }

  for (j = 1; j <=dim; j++)  
    if (rowsol[colsol[j]] != j)
    {
      printf("\n");
      printf("error in col solution j= %d,colsol[j]= %d,rowsol[colsol[j]]= %d\n", j, colsol[j], rowsol[colsol[j]]);
      printf("\n\ndim %5d - press key\n", dim);
      //scanf("%d", &wait);
      break; 
    }

  delete[] matched;
  return;
}


