#include <math.h>

void minusPt(point init, point end, point &result){
	for(int i=0;i<3;i++)
		result[i]=end[i]-init[i];
	}

void addPt(point a, point b, point & result){
	for (int i=0;i<3;i++)
		result[i]=a[i]+b[i];
}

void scalePt(point a, float scale, point &result){
	for(int i=0; i<3; i++)
		result[i]=scale*a[i];
}

void decomposePt(point a, float &scalar, point &direction){
	scalar=0; int i=0;
	for(i=0;i<3;i++)
		scalar+=pow(a[i],2);
	scalar=sqrt(scalar);
	for(i=0;i<3;i++)
		if (scalar!=0)//avoid division by zero -->so not funny
			direction[i]=a[i]/scalar;
		else direction[i]=0;
}

inline 
void sumManyPt(vector <point> a, point &result){
	int i,j;
	for(i=0; i<(int)a.size();i++)
		for(j=0;j<3;j++)
			result[i]+=a[i][j];
}

float dist2Pt(point a, point b){
	int i;
	float result=0;
	for(i=0;i<3;i++)
		result+=pow(a[i]-b[i],2);
	result=sqrt(result);
	return result;
}

int findMinId(vector<float>input){ //find indice of the minimum value from an input array
	int i, id;
	float min=1e20;
	for(i=0;i<(int)input.size();i++)
		if(input[i]<min)
			id=i;
	return id;
}

void badSortId(vector<float> input, vector<int> &idSorted, int num){ //sort in ascending order (min first) and not optimized at all ->run time will be n! until num
	int i,j,k;
	float min=1e20;
	vector <int> flag;
	for(i=0;i<(int)input.size();i++)
		flag.push_back(0);
	for(j=0;j<num;j++){
		for(i=0;i<(int)input.size();i++)
			if(flag[i]==0 && input[i]<=min){
				min=input[i];
				k=i;
			}
			flag[k]=1;
			idSorted.push_back(k);
			min=1e20;
	}
}

//float SphereDist(point a, point b){//calculate spherical distance between 2 points on the same sphere



		

