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
#include <sys/time.h>
#include "vector.h"
#include "vectorAlgebra.h"

typedef struct{
	Vector pos;
	Vector vel;
	Vector fAct;
	Vector acel;
	double masa;
} Body;

double bodiesDist(Body c1, Body c2);

void updF(Body bodies[],int arraySize);

double forceGrav(Body c1, Body c2);

void updAcelOld(Body bodies[], int arraySize);
void updPosOld(Body bodies[], int arraySize, double deltaT);
void updVelOld(Body bodies[], int arraySize, double deltaT);

void updAcelFey(Body bodies[], int arraySize);
void updPosFey(Body bodies[], int arraySize, double deltaT);
void updVelFey(Body bodies[], int arraySize, double deltaT);
void updVelFeyInit(Body bodies[], int arraySize, double deltaT);

double CONSTANTE_GRAVITACION=6.67428e-11;
int DIRECT_METHOD=1;
int FEYNMAN_METHOD=2;


double totalEnergy(Body bodies[], int bodiesSize);
double kinEnergy(Body bodies[], int bodiesSize);
double potEnery(Body bodies[], int bodiesSize);

void printPosData(Body bodies[], int arraySize, FILE *file, int dataType);
void printEnerData(Body bodies[], int arraySize, FILE *file, int dataType);



int main(void){
	/*velocidades y posiciones iniciales de dos particulas*/
    Vector pos1, vel1;	
	Body sol,tie,luna,venus;
	int steps;
	int const bodiesSize=4;
	Body bodies[4];
	double deltaT;
	int cont;
	FILE *file,*enerFile;
	int fileOutInterval;
	int amountOfPoints;	
	int outData;	
    struct timeval timeStart, timeEnd;
    gettimeofday(&timeStart,NULL);
	
    file = fopen("outFey","w"); 
	enerFile = fopen("ener","w");	
	if(file == NULL || enerFile == NULL)  {
        	perror("failed to open some file");
        	return EXIT_FAILURE;
	}

	
	deltaT=16004;
	steps=5200000;
	amountOfPoints=176;	
	fileOutInterval = steps/amountOfPoints;
	
	/*Data Init*/
	/*Sun*/
	pos1.x=0;
	pos1.y=0;
	pos1.z=0;	
	vel1.x=0;
	vel1.y=0;
	vel1.z=0;
	sol.pos=pos1;
	sol.vel=vel1;
	sol.masa=1.989e30;	
	
	/*Earth*/
	pos1.x=1.496e11;
	pos1.y=0;
	pos1.z=0;
	vel1.x=0;
	vel1.y=30000;
	vel1.z=0;
	tie.masa=5.9736e24;	
	tie.pos=pos1;
	tie.vel=vel1;

	/*Moon*/
	pos1.x=tie.pos.x+3.84339e8;
	pos1.y=tie.pos.y;
	pos1.z=tie.pos.z;
	vel1.x=0;
	vel1.y=tie.vel.y+1000;
	vel1.z=0;
	luna.masa=7.3477e22;
	luna.pos=pos1;
	luna.vel=vel1;

	/*venus*/
	pos1.x=1.082e11;
	pos1.y=0;
	pos1.z=0;
	vel1.x=0;
	vel1.y=35021;
	vel1.z=0;
	venus.masa=4.869e24;	
	venus.pos=pos1;
	venus.vel=vel1;
		
	bodies[0]=sol;
	bodies[1]=tie;	
	bodies[2]=luna;
	bodies[3]=venus;
	
	/*init feynman computation*/
	updF(bodies,bodiesSize);
	updVelFeyInit(bodies,bodiesSize,bodiesSize);

	for (cont=0;cont<steps;cont++){
		outData=!(cont%fileOutInterval);
		updPosFey(bodies, bodiesSize, deltaT);
		updF(bodies, bodiesSize);
		updAcelFey(bodies,bodiesSize);
		updVelFey(bodies, bodiesSize, deltaT);
		if(outData){
			printPosData(bodies,bodiesSize,file,FEYNMAN_METHOD);
			printEnerData(bodies,bodiesSize,enerFile,FEYNMAN_METHOD);
		}
		/*
		updAcelOld(bodies,bodiesSize);
		updVelOld(bodies, bodiesSize, deltaT);
		updPosOld(bodies, bodiesSize, deltaT);
		*/
	}


	fclose(file);
	fclose(enerFile);
    gettimeofday(&timeEnd,NULL);
    printf("Total Running Time (secs): %ld.%06ld\n",(timeEnd.tv_usec-timeStart.tv_usec)>0?timeEnd.tv_sec-timeStart.tv_sec:timeEnd.tv_sec-1-timeStart.tv_sec,(timeEnd.tv_usec-timeStart.tv_usec)>0?timeEnd.tv_usec-timeStart.tv_usec:timeStart.tv_usec-timeEnd.tv_usec);
	return 0;
}

double bodiesDist(Body c1, Body c2){
	return vectoresDist(c1.pos,c2.pos);
}

double forceGrav(Body c1, Body c2){
	double toReturn;
	toReturn= c1.masa*c2.masa*CONSTANTE_GRAVITACION/pow(bodiesDist(c1,c2),2);
	return toReturn;
}

void updF(Body bodies[], int arraySize){
	int j,k;
	double force;
	Vector dirR;
	double distancia;

	/*
	Se limpian las forces de la iteracion anterior
	*/
	for(j=0;j<arraySize;j++){
		bodies[j].fAct.x=0;
		bodies[j].fAct.y=0;
		bodies[j].fAct.z=0;
	}	

	for(j=0;j<arraySize;j++){
		for (k=j+1;k<arraySize;k++){
			force = forceGrav(bodies[j],bodies[k]);
			dirR=vectoresResta(bodies[j].pos,bodies[k].pos);
			distancia=vectorModulo(dirR);	
			
			/*es un menos por que el vector apunta de k hacia j*/
			bodies[j].fAct.x+=-force*dirR.x/distancia;	
			bodies[j].fAct.y+=-force*dirR.y/distancia;	
			bodies[j].fAct.z+=-force*dirR.z/distancia;	
	
			bodies[k].fAct.x+=force*dirR.x/distancia;	
			bodies[k].fAct.y+=force*dirR.y/distancia;	
			bodies[k].fAct.z+=force*dirR.z/distancia;	
		}	
	}
}

void updAcelOld(Body bodies[], int arraySize){
	int i;
	for(i = 0; i<arraySize;i++){
		bodies[i].acel.x=bodies[i].fAct.x/bodies[i].masa;
		bodies[i].acel.y=bodies[i].fAct.y/bodies[i].masa;
		bodies[i].acel.z=bodies[i].fAct.z/bodies[i].masa;
	}
}

void updAcelFey(Body bodies[], int arraySize){
	int i;
	for(i = 0; i<arraySize;i++){
		bodies[i].acel.x=bodies[i].fAct.x/bodies[i].masa;
		bodies[i].acel.y=bodies[i].fAct.y/bodies[i].masa;
		bodies[i].acel.z=bodies[i].fAct.z/bodies[i].masa;
	}
}

void updVelFey(Body bodies[], int arraySize, double deltaT){
	int i;
	Body c;
	for (i=0;i<arraySize;i++){
		c = bodies[i];
		bodies[i].vel.x=bodies[i].vel.x+bodies[i].acel.x*deltaT;
		bodies[i].vel.y=bodies[i].vel.y+bodies[i].acel.y*deltaT;
		bodies[i].vel.z=bodies[i].vel.z+bodies[i].acel.z*deltaT;
	}
}

void updVelFeyInit(Body bodies[], int arraySize, double deltaT){
	int i;
	Body c;
	for (i=0;i<arraySize;i++){
		c = bodies[i];
		bodies[i].vel.x=bodies[i].vel.x+bodies[i].acel.x*deltaT/2;
		bodies[i].vel.y=bodies[i].vel.y+bodies[i].acel.y*deltaT/2;
		bodies[i].vel.z=bodies[i].vel.z+bodies[i].acel.z*deltaT/2;
	}
}

void updPosFey(Body bodies[], int arraySize, double deltaT){
	int i;
	Body c;
	for (i=0;i<arraySize;i++){
		c = bodies[i];
	
		bodies[i].pos.x+=bodies[i].vel.x*deltaT;
		bodies[i].pos.y+=bodies[i].vel.y*deltaT;
		bodies[i].pos.z+=bodies[i].vel.z*deltaT;
		
	}
}

void printPosData(Body bodies[], int arraySize, FILE *file, int dataType){
	char out[1000]="",buffer[60];
	int i = 0;	
	for (;i<arraySize;i++){
		sprintf(buffer,"%11.3e \t %11.3e \t",bodies[i].pos.x,bodies[i].pos.y);
		strcat(out,buffer);
	}
	strcat(out,"\n");
	fwrite(out,1,strlen(out),file);

}

void printEnerData(Body bodies[], int arraySize, FILE *file, int dataType){
	char out[1000]="",buffer[60];
	sprintf(buffer,"%11.3e \t ",totalEnergy(bodies,arraySize));
	strcat(out,buffer);
	strcat(out,"\n");
	fwrite(out,1,strlen(out),file);

}

void updPosOld(Body bodies[], int arraySize, double deltaT){
	int i;
	Body c;
	char out[1000]="",buffer[50];	
	
	for (i=0;i<arraySize;i++){
		c = bodies[i];
	
		bodies[i].pos.x+=bodies[i].vel.x*deltaT+bodies[i].acel.x*deltaT*deltaT/2;
		bodies[i].pos.y+=bodies[i].vel.y*deltaT+bodies[i].acel.y*deltaT*deltaT/2;
		bodies[i].pos.z+=bodies[i].vel.z*deltaT+bodies[i].acel.z*deltaT*deltaT/2;
		sprintf(buffer,"%11.3e \t %11.3e \t",bodies[i].pos.x,bodies[i].pos.y);
		strcat(out,buffer);
	}
	printf("%s  \n",out);
}


void updVelOld(Body bodies[], int arraySize, double deltaT){
	int i;
	Body c;
	for (i=0;i<arraySize;i++){
		c = bodies[i];
		bodies[i].vel.x+=bodies[i].acel.x*deltaT;
		bodies[i].vel.y+=bodies[i].acel.y*deltaT;
		bodies[i].vel.z+=bodies[i].acel.z*deltaT;
	}
}

double totalEnergy(Body bodies[], int bodiesSize){
	return kinEnergy(bodies,bodiesSize)+potEnery(bodies,bodiesSize);	
}



double kinEnergy(Body bodies[], int bodiesSize){
	int i=0;
	double ener=0;
	Body cuerpo;
	double vel=0;
	for(;i<bodiesSize;i++){
		cuerpo=bodies[i];
		vel=vectorModulo(cuerpo.vel);
		ener+=cuerpo.masa*vel*vel/2;
	}
	return ener;
}

double potEnery(Body bodies[], int bodiesSize){
	int i=0,j=0;
	double ener=0;
	double dist;
	double masa;
	for(;i<bodiesSize;i++){
		masa=bodies[i].masa;
		for(j=i+1;j<bodiesSize;j++){
			dist=vectoresDist(bodies[i].pos,bodies[j].pos);
			ener+=-CONSTANTE_GRAVITACION*masa*bodies[j].masa/dist;
		}
	}
	return ener;
}

Vector totalAngMom(Body bodies[], int bodiesSize){
	Vector toRet;
	toRet.x=0;
	toRet.y=0;
	toRet.z=0;
		
	return toRet;	
}

