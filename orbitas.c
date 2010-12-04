/*
Programa para calcular orbitas de planetas:
TODO LIST
_ Mejorar los calculos de forces y/o velocidades (Feynman tiene una linda explicacion)
_ Hacer que se lean de ficheros externos los datos de los problemas
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

typedef struct{
	double xEnd;
	double yEnd;
	double zEnd;
	double xStart;
	double yStart;
	double zStart;
} Cuadrante;

typedef struct{
	Vector pos;
	Vector vel;
	Vector fAct;
	Vector acel;
	double masa;
    Cuadrante *cuadrante;
} Body;

Cuadrante space(Body bodies[], int qBodies);

Cuadrante space(Body bodies[], int qBodies){
    Cuadrante c;
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

    return c;
}


/*Slices representa en cuantos trozos se partira cada dimension a la hora de hacer los cuadrantesl*/
void cuadrantePartition(int slices, Cuadrante c[], Body bodies[], int qBodies); 

void cuadrantePartition(int slices, Cuadrante c[], Body bodies[], int qBodies){
   Cuadrante spaceT;
   spaceT = space(bodies, qBodies);

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

	deltaT = 540;
	steps = 420000;
	amountOfPoints = 77037;	
	fileOutInterval = steps/amountOfPoints;

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
		sprintf(buffer, "%11.3e \t %11.3e \t", bodies[i].pos.x/1000000, bodies[i].pos.y/1000000);
		strcat(out, buffer);
	}
	strcat(out, "\n");
	fwrite(out, 1, strlen(out), file);

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



