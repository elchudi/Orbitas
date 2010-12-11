/*
Programa para calcular orbitas de planetas:
TODO LIST
_ Hacer que se resuelvan mediante autovectores y/o sistemas de equaciones
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vector.h"
#include "vectorAlgebra.h"

#define qMax 300

typedef struct tOctant{
	double xEnd;
	double yEnd;
	double zEnd;
	double xStart;
	double yStart;
	double zStart;
    struct tOctant *parent;
    struct tOctant *childs[8];
    int isLeaf;
    int isRoot;
} Octant  ;

typedef struct tBody{
	Vector pos;
	Vector vel;
	Vector fAct;
	Vector acel;
	double masa;
} Body;

/*Octant constructor !!! is this the way to do it? (i.e., if I need to create and array of poiters of structs in a function, should i do this or is ther such a thing as 'new struct' */
struct tOctant* getOctant (double xS, double yS, double zS, double xE, double yE, double zE);
struct tOctant* getOctant (double xS, double yS, double zS, double xE, double yE, double zE){
    Octant *toRet;
    toRet = malloc(sizeof(Octant));
    (*toRet).xStart = xS;
    (*toRet).yStart = yS;
    (*toRet).zStart = zS;
    (*toRet).xEnd = xE;
    (*toRet).yEnd = yE,
    (*toRet).zEnd = zE;
   
    return toRet;
}

Octant space(Body bodies[], int qBodies);

Octant space(Body bodies[], int qBodies){
    Octant c;
    int i;
    c.xStart = bodies[0].pos.x;
    c.yStart = bodies[0].pos.y;
    c.zStart = bodies[0].pos.z;
    c.xEnd = bodies[0].pos.x;
    c.yEnd = bodies[0].pos.y;
    c.zEnd = bodies[0].pos.z;


    for(i = 1; i < qBodies; i++){
        if( c.xStart > bodies[i].pos.x )   c.xStart = bodies[i].pos.x;
        if( c.yStart > bodies[i].pos.y )   c.yStart = bodies[i].pos.y;
        if( c.zStart > bodies[i].pos.z )   c.zStart = bodies[i].pos.z;
        if( c.xEnd < bodies[i].pos.x )   c.xEnd = bodies[i].pos.x;
        if( c.yEnd < bodies[i].pos.y )   c.yEnd = bodies[i].pos.y;
        if( c.zEnd < bodies[i].pos.z )   c.zEnd = bodies[i].pos.z;
        
    }
    /*Some breathing air for the space*/
    if(c.xEnd > 0){ c.xEnd = c.xEnd*1.1;}else{ c.xEnd = c.xEnd*0.9;}
    if(c.yEnd > 0){ c.yEnd = c.yEnd*1.1;}else{ c.yEnd = c.yEnd*0.9;}
    if(c.zEnd > 0){ c.zEnd = c.zEnd*1.1;}else{ c.zEnd = c.zEnd*0.9;}
    if(c.xStart < 0){ c.xStart = c.xStart*1.1;}else{ c.xStart = c.xStart*0.9;}
    if(c.yStart < 0){ c.yStart = c.yStart*1.1;}else{ c.yStart = c.yStart*0.9;}
    if(c.zStart < 0){ c.zStart = c.zStart*1.1;}else{ c.zStart = c.zStart*0.9;}

    return c;
}


void octantPartition(Body bodies[], int qBodies); 

void octantPartition(Body bodies[], int qBodies){
    Octant spaceT;
    int i = 0, j = 0, k = 0, slices = 0, cont = 0;
    double xLen = 0, yLen = 0, zLen = 0;
    /*Total space to partition*/
    spaceT = space(bodies, qBodies);
    /*Size of each octant*/
    slices = 2;
    if(slices != 0){
        xLen = (spaceT.xEnd - spaceT.xStart)/slices;
        yLen = (spaceT.yEnd - spaceT.yStart)/slices;
        zLen = (spaceT.zEnd - spaceT.zStart)/slices;
    }else{
        xLen = spaceT.xEnd - spaceT.xStart;
        yLen = spaceT.yEnd - spaceT.yStart;
        zLen = spaceT.zEnd - spaceT.zStart;
    }
    cont = 0;
    for(i = 0; i < slices; i++){
        for (j = 0; j < slices; j++){
            for (k = 0; k < slices; k++){
                
                /*
                aux.xStart = spaceT.xStart + xLen*i;    
                aux.yStart = spaceT.yStart + yLen*j;    
                aux.zStart = spaceT.zStart + zLen*k;    
                aux.xEnd = spaceT.xStart + xLen*(i+1);
                aux.yEnd = spaceT.yStart + yLen*(j+1);
                aux.zEnd = spaceT.zStart + zLen*(k+1);
                */
                spaceT.childs[cont] = getOctant(spaceT.xStart + xLen*i, spaceT.yStart + yLen*j, spaceT.zStart + zLen*k, spaceT.xStart + xLen*(i+1), spaceT.yStart + yLen*(j+1), spaceT.zStart + zLen*(k+1));
                cont++;
            }
        }
    } 

} 



double bodiesDist(Body c1, Body c2);

void updF(Body bodies[],int qBodies);

double forceGrav(Body c1, Body c2);

void updAcelOld(Body bodies[], int qBodies);
void updPosOld(Body bodies[], int qBodies, double deltaT);
void updVelOld(Body bodies[], int qBodies, double deltaT);

void updAcelFey(Body bodies[], int qBodies);
void updPosFey(Body bodies[], int qBodies, double deltaT);
void updVelFey(Body bodies[], int qBodies, double deltaT);
void updVelFeyInit(Body bodies[], int qBodies, double deltaT);

double CONSTANTE_GRAVITACION=6.67428e-11;
int DIRECT_METHOD=1;
int FEYNMAN_METHOD=2;

void loadExtData(Body bodies[],int *qBodies, FILE *file);

double totalEnergy(Body bodies[], int qBodies);
double kinEnergy(Body bodies[], int qBodies);
double potEnery(Body bodies[], int qBodies);

void printPosData(Body bodies[], int qBodies, FILE *file, int dataType);
void printEnerData(Body bodies[], int qBodies, FILE *file, int dataType);
void printPlotingFile(int qBodies );



int main(void){
	int steps;
	int qBodies = 0;
	Body bodies[qMax];
	double deltaT;
	int cont;
	FILE *file, *enerFile, *fileData;
	int fileOutInterval;
	int amountOfPoints;	
    clock_t start, end;
    double runTime;
    start = clock();
     
    file = fopen("outFey", "w"); 
	enerFile = fopen("ener", "w");	
	if(file == NULL || enerFile == NULL)  {
        	perror("failed to open some file");
        	return EXIT_FAILURE;
	}
		
	fileData = fopen("initData", "r");
    loadExtData(bodies, &qBodies, fileData);

	deltaT = 4400;
	steps = 220000;
	amountOfPoints = 37037;	
	fileOutInterval = steps/amountOfPoints;
    printPlotingFile(qBodies);
	/*init feynman computation*/
	updF(bodies, qBodies);
	updVelFeyInit(bodies, qBodies, deltaT);

	for (cont=0;cont<steps;cont++){
		updPosFey(bodies, qBodies, deltaT);
		updF(bodies, qBodies);
		updAcelFey(bodies, qBodies);
		updVelFey(bodies, qBodies, deltaT);
		if(!(cont%fileOutInterval)){
			printPosData(bodies, qBodies, file, FEYNMAN_METHOD);
			printEnerData(bodies, qBodies, enerFile, FEYNMAN_METHOD);
		}
		/*
		updAcelOld(bodies, qBodies);
		updVelOld(bodies, qBodies, deltaT);
		updPosOld(bodies, qBodies, deltaT);
		*/
	}


	fclose(file);
	fclose(enerFile);
	fclose(fileData);
	end = clock();
    runTime = (end - start) / (double) CLOCKS_PER_SEC ;
    printf ("Run time is %g seconds\n", runTime);
    
    return 0;
}

double bodiesDist(Body c1, Body c2){
	return vectoresDist(&c1.pos, &c2.pos);
}

double forceGrav(Body c1, Body c2){
	double toReturn;
	toReturn = c1.masa*c2.masa*CONSTANTE_GRAVITACION/pow(bodiesDist(c1, c2), 2);
	return toReturn;
}

void updF(Body bodies[], int qBodies){
	int j,k;
	double force;
	Vector dirR;
	double distancia;

	/*
	Se limpian las forces de la iteracion anterior
	*/
	for(j = 0; j < qBodies; j++){
		bodies[j].fAct.x = 0;
		bodies[j].fAct.y = 0;
		bodies[j].fAct.z = 0;
	}	
	for(j=0;j< qBodies;j++){
		for (k = j + 1; k < qBodies; k++){
			force = forceGrav(bodies[j],bodies[k]);
    		dirR = vectoresResta(&bodies[j].pos,&bodies[k].pos);
			distancia = vectorModulo(&dirR);	
			
			/*es un menos por que el vector apunta de k hacia j*/
			bodies[j].fAct.x += -force*dirR.x/distancia;	
			bodies[j].fAct.y += -force*dirR.y/distancia;	
			bodies[j].fAct.z += -force*dirR.z/distancia;	
	
			bodies[k].fAct.x += force*dirR.x/distancia;	
			bodies[k].fAct.y += force*dirR.y/distancia;	
			bodies[k].fAct.z += force*dirR.z/distancia;	
		}	
	}
    
}

void updAcelOld(Body bodies[], int qBodies){
	int i;
	for(i = 0; i < qBodies; i++){
		bodies[i].acel.x = bodies[i].fAct.x/bodies[i].masa;
		bodies[i].acel.y = bodies[i].fAct.y/bodies[i].masa;
		bodies[i].acel.z = bodies[i].fAct.z/bodies[i].masa;
	}
}

void updAcelFey(Body bodies[], int qBodies){
	int i;
	for(i = 0; i < qBodies; i++){
		bodies[i].acel.x = bodies[i].fAct.x/bodies[i].masa;
		bodies[i].acel.y = bodies[i].fAct.y/bodies[i].masa;
		bodies[i].acel.z = bodies[i].fAct.z/bodies[i].masa;
	}
}

void updVelFey(Body bodies[], int qBodies, double deltaT){
	int i;
	Body c;
	for (i = 0; i < qBodies; i++){
		c = bodies[i];
		bodies[i].vel.x = bodies[i].vel.x + bodies[i].acel.x*deltaT;
		bodies[i].vel.y = bodies[i].vel.y + bodies[i].acel.y*deltaT;
		bodies[i].vel.z = bodies[i].vel.z + bodies[i].acel.z*deltaT;
	}
}

void updVelFeyInit(Body bodies[], int qBodies, double deltaT){
	int i;
	Body c;
	for (i = 0; i < qBodies; i++){
		c = bodies[i];
		bodies[i].vel.x = bodies[i].vel.x + bodies[i].acel.x*deltaT/2;
		bodies[i].vel.y = bodies[i].vel.y + bodies[i].acel.y*deltaT/2;
		bodies[i].vel.z = bodies[i].vel.z + bodies[i].acel.z*deltaT/2;
	}
}

void updPosFey(Body bodies[], int qBodies, double deltaT){
	int i;
	Body c;
	for (i = 0; i < qBodies; i++){
		c = bodies[i];
	
		bodies[i].pos.x += bodies[i].vel.x*deltaT;
		bodies[i].pos.y += bodies[i].vel.y*deltaT;
		bodies[i].pos.z += bodies[i].vel.z*deltaT;
		
	}
}

void printPosData(Body bodies[], int qBodies, FILE *file, int dataType){
	char out[1000] = "", buffer[60];
	int i = 0;	
	for (; i < qBodies; i++){

        /*Sin centro de referencia
		sprintf(buffer, "%11.3e \t %11.3e \t", bodies[i].pos.x/1000000, bodies[i].pos.y/1000000);
        */
        /*visto desde la tierra
		sprintf(buffer, "%11.3e \t %11.3e \t", bodies[i].pos.x/1000000-bodies[3].pos.x/1000000, bodies[i].pos.y/1000000-bodies[3].pos.y/1000000);
        */
		/*visto desde el sol*/
        sprintf(buffer, "%11.3e \t %11.3e \t", bodies[i].pos.x/1000000-bodies[0].pos.x/1000000, bodies[i].pos.y/1000000-bodies[0].pos.y/1000000);
	    strcat(out, buffer);
	}
	strcat(out, "\n");
	fwrite(out, 1, strlen(out), file);

}

void printPlotingFile(int qBodies ){
	char out[1000] = "", buffer[60];
    int i;
    FILE *file;
    file = fopen("orbitas.p", "w"); 
	strcat(out, "plot 'outFey' using 1:2 with lines");
	for (i=2; i < qBodies; i++){
		sprintf(buffer, ", 'outFey' using %d:%d with lines",i*2-1,2*i);
		strcat(out, buffer);
	}
	strcat(out, "\n");
	fwrite(out, 1, strlen(out), file);
    fclose(file);
}

void printEnerData(Body bodies[], int qBodies, FILE *file, int dataType){
	char out[1000] = "", buffer[60];
	sprintf(buffer,"%11.3e \t ", totalEnergy(bodies, qBodies));
	strcat(out, buffer);
	strcat(out, "\n");
	fwrite(out, 1, strlen(out), file);

}

void updPosOld(Body bodies[], int qBodies, double deltaT){
	int i;
	Body c;
	char out[1000] = "", buffer[50];	
	
	for (i=0;i< qBodies;i++){
		c = bodies[i];
	
		bodies[i].pos.x += bodies[i].vel.x*deltaT + bodies[i].acel.x*deltaT*deltaT/2;
		bodies[i].pos.y += bodies[i].vel.y*deltaT + bodies[i].acel.y*deltaT*deltaT/2;
		bodies[i].pos.z += bodies[i].vel.z*deltaT + bodies[i].acel.z*deltaT*deltaT/2;
		sprintf(buffer,"%11.3e \t %11.3e \t", bodies[i].pos.x, bodies[i].pos.y);
		strcat(out,buffer);
	}
	printf("%s  \n",out);
}

void updVelOld(Body bodies[], int qBodies, double deltaT){
	int i;
	Body c;
	for (i = 0; i < qBodies; i++){
		c = bodies[i];
		bodies[i].vel.x += bodies[i].acel.x*deltaT;
		bodies[i].vel.y += bodies[i].acel.y*deltaT;
		bodies[i].vel.z += bodies[i].acel.z*deltaT;
	}
}

double totalEnergy(Body bodies[], int qBodies){
	return kinEnergy(bodies, qBodies) + potEnery(bodies, qBodies);	
}

double kinEnergy(Body bodies[], int qBodies){
	int i = 0;
	double ener = 0;
	Body cuerpo;
	double vel = 0;
	for(;i < qBodies; i++){
		cuerpo = bodies[i];
		vel = vectorModulo(&cuerpo.vel);
		ener += cuerpo.masa*vel*vel/2;
	}
	return ener;
}

double potEnery(Body bodies[], int qBodies){
	int i = 0,j = 0;
	double ener = 0;
	double dist;
	double masa;
	for(;i < qBodies; i++){
		masa = bodies[i].masa;
		for(j = i + 1;j < qBodies; j++){
			dist = vectoresDist(&bodies[i].pos,&bodies[j].pos);
			ener += -CONSTANTE_GRAVITACION*masa*bodies[j].masa/dist;
		}
	}
	return ener;
}

Vector totalAngMom(Body bodies[], int qBodies){
	Vector toRet;
	toRet.x = 0;
	toRet.y = 0;
	toRet.z = 0;
		
	return toRet;	
}

void loadExtData(Body bodies[], int *qBodies, FILE *file){
    int i = 0;
    rewind(file);
    fscanf(file, "%d\n", qBodies);
    for (;i < *qBodies; i++){
        fscanf(file, "%lf%lf%lf%lf%lf%lf%lf\n", &bodies[i].masa, &bodies[i].pos.x, &bodies[i].pos.y, &bodies[i].pos.z, &bodies[i].vel.x, &bodies[i].vel.y, &bodies[i].vel.z);
        /*printf("%g %g %g %g %g %g %g\n",bodies[i].masa,bodies[i].pos.x,bodies[i].pos.y,bodies[i].pos.z,bodies[i].vel.x,bodies[i].vel.y,bodies[i].vel.z);*/
    }
}
