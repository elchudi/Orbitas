#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "vectorAlgebra.h"



double vectorCoorModulo(double x1, double y1, double z1){
	float dist;
	dist=0.f;
	dist = sqrt(x1*x1+y1*y1+z1*z1);
	return dist;	
}

double vectorModulo (Vector vec){
	return vectorCoorModulo(vec.x,vec.y,vec.z);
}

Vector vectoresResta(Vector vec1, Vector vec2){
	Vector resul;
	resul.x=vec1.x-vec2.x;
	resul.y=vec1.y-vec2.y;
	resul.z=vec1.z-vec2.z;
	return resul;
}

double vectoresDist(Vector vec1, Vector vec2){
	return vectorModulo(vectoresResta(vec1,vec2));
}

Vector crossProd(Vector x, Vector y){
	Vector toRet;
	toRet.x=x.y*y.z-x.z*y.y;
	toRet.y=x.z*y.x-x.x*y.z;
	toRet.z=x.x*y.y-x.y*y.x;	
	return toRet;
}

Vector vectorSum(Vector x, Vector y){
	Vector toRet;
	toRet.x=x.x+y.x;
	toRet.y=x.y+y.y;
	toRet.z=x.z+y.z;
	return toRet;
}


/*esto como mierda se usa eh?*/
void vectorInitToZero(Vector *x){
	(*x).x=0;
	(*x).y=0;
	(*x).z=0;
}

