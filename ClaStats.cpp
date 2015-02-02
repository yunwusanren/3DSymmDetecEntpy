
/** ClaStats parses a .cla file, verifies that it is correctly formated, and then prints
    out statistics about the classification.

    Usage
    ./CreateClaOverview file.cla

    Copyright 2003 Princeton University

*/
#include <stdafx.h>
#include <stdio.h>

#include "PSBClaParse.h"


extern PSBCategoryList* ClassFileParse(char *classfilename)
{
 PSBCategoryList* categories;
  PSBCategory* category;
  int totalClasses, nonEmpty, i, totalModels, largestClass;

  categories = parseFile(classfilename);

  totalClasses = categories->_numCategories;
  largestClass = totalModels = nonEmpty = 0;

  for(i = 0; i < categories->_numCategories; ++i){
    category = categories->_categories[i];
    if (category->_numModels > 0){
      nonEmpty++;
      totalModels+= category->_numModels; 
      if (category->_numModels > largestClass){
        largestClass = category->_numModels;
      }
    }
  }
/*
  printf("Classification file: %s\n",classfilename);
  printf("Num classes: %d\n", totalClasses);
  printf("Num non-empty classes: %d\n", nonEmpty);
  printf("Num models: %d\n", totalModels); 
  printf("Average class size (non-empty): %g\n", totalModels/(float)nonEmpty);
  printf("Largest class size: %d\n", largestClass);
*/
  return categories;
}
