#include <stdio.h>
#include <math.h>
#include "vector.h"
#include "vectorAlgebra.h"

double vectorCoorModulo(double *x1, double *y1, double *z1){
	return sqrt((*x1)*(*x1)+(*y1)*(*y1)+(*z1)*(*z1));
}

double vectorModulo (Vector *vec){
	return vectorCoorModulo(&(*vec).x,&(*vec).y,&(*vec).z);
}

Vector vectoresResta(Vector *vec1, Vector *vec2){
	Vector resul;
	resul.x = (*vec1).x - (*vec2).x;
	resul.y = (*vec1).y - (*vec2).y;
	resul.z = (*vec1).z - (*vec2).z;
	return resul;
}

/*se puede d*/
double vectoresDist(Vector *vec1, Vector *vec2){
    Vector aux =vectoresResta(vec1,vec2);
    return vectorModulo(&aux);
}

Vector crossProd(Vector *x, Vector *y){
	Vector toRet;
	toRet.x = (*x).y*(*y).z - (*x).z*(*y).y;
	toRet.y = (*x).z*(*y).x - (*x).x*(*y).z;
	toRet.z = (*x).x*(*y).y - (*x).y*(*y).x;	
	return toRet;
}

Vector vectorSum(Vector *x, Vector *y){
	Vector toRet;
	toRet.x = (*x).x + (*y).x;
	toRet.y = (*x).y + (*y).y;
	toRet.z = (*x).z + (*y).z;
	return toRet;
}

void vectorInitToZero(Vector *x){
	(*x).x = 0;
	(*x).y = 0;
	(*x).z = 0;
}

