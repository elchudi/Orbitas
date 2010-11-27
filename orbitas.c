/*
Programa para calcular orbitas de planetas:
TODO LIST
_ Mejorar los calculos de fuerzas y/o velocidades (Feynman tiene una linda explicacion)
_ Hacer que se lean de ficheros externos los datos de los problemas
_ Hacer que se resuelvan mediante autovectores y/o sistemas de equaciones


*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

typedef struct{
	double x;
	double y;
	double z;
} Vector; 

typedef struct{
	Vector pos;
	Vector vel;
	Vector fAct;
	Vector acel;
	double masa;
} Cuerpo;

Vector vectoresResta(Vector vec1, Vector vec2);

double cuerposDist(Cuerpo c1, Cuerpo c2);
double vectoresDist(Vector vec1, Vector vec2);
double vectorModulo(Vector vec);

void holaMundo(void);

double vectorCoorModulo(double x1, double y1, double z1);

void updF(Cuerpo cuerpos[],int arraySize);


double fuerzaGrav(Cuerpo c1, Cuerpo c2);



void updAcelOld(Cuerpo cuerpos[], int arraySize);
void updPosOld(Cuerpo cuerpos[], int arraySize, double deltaT);
void updVelOld(Cuerpo cuerpos[], int arraySize, double deltaT);

void updAcelFey(Cuerpo cuerpos[], int arraySize);
void updPosFey(Cuerpo cuerpos[], int arraySize, double deltaT, FILE *file ,int print);
void updVelFey(Cuerpo cuerpos[], int arraySize, double deltaT);
void updVelFeyInit(Cuerpo cuerpos[], int arraySize, double deltaT);

double CONSTANTE_GRAVITACION=6.67428e-11;



int main(void){
	/*velocidades y posiciones iniciales de dos particulas*/
	Vector pos1, vel1;	
	Cuerpo sol,tie,luna,venus;
	Cuerpo cuerpos[3];
	int steps,cuerposSize;
	double deltaT;
	int cont;
	FILE *file;
	int fileOutInterval;
	int amountOfPoints;	

	file = fopen("outFey","w"); 
	if(file == NULL) {
        	perror("failed to open sample.txt");
        	return EXIT_FAILURE;
	}

	
	deltaT=16004;
	steps=1200000;
	amountOfPoints=30000;	
	fileOutInterval = steps/amountOfPoints;
	


	

	/*inicializar los datos*/
	/*cuerpo 1*/
	pos1.x=0;
	pos1.y=0;
	pos1.z=0;	
	vel1.x=0;
	vel1.y=0;
	vel1.z=0;
	sol.pos=pos1;
	sol.vel=vel1;
	sol.masa=1.989e30;	
	

	/*cuerpo 2*/
	pos1.x=1.496e11;
	pos1.y=0;
	pos1.z=0;
	vel1.x=0;
	vel1.y=30000;
	vel1.z=0;
	tie.masa=5.9736e24;	
	tie.pos=pos1;
	tie.vel=vel1;
		

	/*luna*/
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
		
	

	cuerpos[0]=sol;
	cuerpos[1]=tie;	
	cuerpos[2]=luna;
	cuerpos[3]=venus;
	cuerposSize=4;	
	
	/*
	printf("la distancia entre los cuerpos es %g \n",cuerposDist(sol,tie));
	printf("la fuerza entre los cuerpos es %g \n",fuerzaGrav(sol,tie));
	printf("la distancia entre vectores es %g \n",vectoresDist(pos1,pos2));
	printf("objeto con direccion %f",*cuerpos);
	printf("masa de 1 %f\n",cuerpos[0].masa);
	printf("masa de 1 %f\n",cuerpos[0].masa);
	*/

	/*init feynman computation*/
	updF(cuerpos,cuerposSize);
	updVelFeyInit(cuerpos,cuerposSize,cuerposSize);

	for (cont=0;cont<steps;cont++){
		updPosFey(cuerpos, cuerposSize, deltaT, file,!(cont%fileOutInterval));
		updF(cuerpos, cuerposSize);
		updAcelFey(cuerpos,cuerposSize);
		updVelFey(cuerpos, cuerposSize, deltaT);

		/*
		updAcelOld(cuerpos,cuerposSize);
		updVelOld(cuerpos, cuerposSize, deltaT);
		updPosOld(cuerpos, cuerposSize, deltaT);
		*/
	}


	fclose(file);
	return 0;
}




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

double cuerposDist(Cuerpo c1, Cuerpo c2){
	return vectoresDist(c1.pos,c2.pos);
}

double fuerzaGrav(Cuerpo c1, Cuerpo c2){
	double toReturn;
	toReturn= c1.masa*c2.masa*CONSTANTE_GRAVITACION/pow(cuerposDist(c1,c2),2);
	return toReturn;
}

void updF(Cuerpo cuerpos[], int arraySize){
	int j,k;
	double force;
	Vector dirR;
	double distancia;

	/*	
	printf("objeto con direccion %f",*cuerpos);
	arraySize=sizeof(cuerpos)/sizeof(cuerpos[0]);
	printf("sizeOf(cuerpos) %d\n",sizeof(cuerpos));
	printf("sizeof(cuerpos[0]) %d\n",sizeof(cuerpos[0]));
	printf("sizeof(*cuerpos) %d\n",sizeof(*cuerpos));
	printf("arraySize %d\n",arraySize);
	*/
	/*
	Se limpian las fuerzas de la iteracion anterior
	*/
	for(j=0;j<arraySize;j++){
		cuerpos[j].fAct.x=0;
		cuerpos[j].fAct.y=0;
		cuerpos[j].fAct.z=0;
	}	

	for(j=0;j<arraySize;j++){
		for (k=j+1;k<arraySize;k++){
			force = fuerzaGrav(cuerpos[j],cuerpos[k]);
			dirR=vectoresResta(cuerpos[j].pos,cuerpos[k].pos);
			distancia=vectorModulo(dirR);	
			
			/*es un menos por que el vector apunta de k hacia j*/			
			cuerpos[j].fAct.x+=-force*dirR.x/distancia;	
			cuerpos[j].fAct.y+=-force*dirR.y/distancia;	
			cuerpos[j].fAct.z+=-force*dirR.z/distancia;	
	
			cuerpos[k].fAct.x+=force*dirR.x/distancia;	
			cuerpos[k].fAct.y+=force*dirR.y/distancia;	
			cuerpos[k].fAct.z+=force*dirR.z/distancia;	
		}	
	}
}


void updAcelOld(Cuerpo cuerpos[], int arraySize){
	int i;
	
	for(i = 0; i<arraySize;i++){
		cuerpos[i].acel.x=cuerpos[i].fAct.x/cuerpos[i].masa;
		cuerpos[i].acel.y=cuerpos[i].fAct.y/cuerpos[i].masa;
		cuerpos[i].acel.z=cuerpos[i].fAct.z/cuerpos[i].masa;
	}
}

void updAcelFey(Cuerpo cuerpos[], int arraySize){
	int i;
	
	for(i = 0; i<arraySize;i++){
		cuerpos[i].acel.x=cuerpos[i].fAct.x/cuerpos[i].masa;
		cuerpos[i].acel.y=cuerpos[i].fAct.y/cuerpos[i].masa;
		cuerpos[i].acel.z=cuerpos[i].fAct.z/cuerpos[i].masa;
	}
}

void updVelFey(Cuerpo cuerpos[], int arraySize, double deltaT){
	int i;
	Cuerpo c;
	for (i=0;i<arraySize;i++){
		c = cuerpos[i];
		cuerpos[i].vel.x=cuerpos[i].vel.x+cuerpos[i].acel.x*deltaT;
		cuerpos[i].vel.y=cuerpos[i].vel.y+cuerpos[i].acel.y*deltaT;
		cuerpos[i].vel.z=cuerpos[i].vel.z+cuerpos[i].acel.z*deltaT;
	}
}

void updVelFeyInit(Cuerpo cuerpos[], int arraySize, double deltaT){
	int i;
	Cuerpo c;
	for (i=0;i<arraySize;i++){
		c = cuerpos[i];
		cuerpos[i].vel.x=cuerpos[i].vel.x+cuerpos[i].acel.x*deltaT/2;
		cuerpos[i].vel.y=cuerpos[i].vel.y+cuerpos[i].acel.y*deltaT/2;
		cuerpos[i].vel.z=cuerpos[i].vel.z+cuerpos[i].acel.z*deltaT/2;
	}
}

void updPosFey(Cuerpo cuerpos[], int arraySize, double deltaT, FILE *file, int print){
	int i;
	Cuerpo c;
	char out[1000]="",buffer[50];	
	
	for (i=0;i<arraySize;i++){
		c = cuerpos[i];
	
		cuerpos[i].pos.x+=cuerpos[i].vel.x*deltaT;
		cuerpos[i].pos.y+=cuerpos[i].vel.y*deltaT;
		cuerpos[i].pos.z+=cuerpos[i].vel.z*deltaT;
		if(print){
			sprintf(buffer,"%11.3e \t %11.3e \t",cuerpos[i].pos.x,cuerpos[i].pos.y);
			strcat(out,buffer);
		}
		
	}
	if(print){
		strcat(out,"\n");
		fwrite(out,1,strlen(out),file);
	}
 	/*	printf("%s  \n",out);*/
}

void updPosOld(Cuerpo cuerpos[], int arraySize, double deltaT){
	int i;
	Cuerpo c;
	char out[1000]="",buffer[50];	
	
	for (i=0;i<arraySize;i++){
		c = cuerpos[i];
	
		cuerpos[i].pos.x+=cuerpos[i].vel.x*deltaT+cuerpos[i].acel.x*deltaT*deltaT/2;
		cuerpos[i].pos.y+=cuerpos[i].vel.y*deltaT+cuerpos[i].acel.y*deltaT*deltaT/2;
		cuerpos[i].pos.z+=cuerpos[i].vel.z*deltaT+cuerpos[i].acel.z*deltaT*deltaT/2;
		sprintf(buffer,"%11.3e \t %11.3e \t",cuerpos[i].pos.x,cuerpos[i].pos.y);
		strcat(out,buffer);
	}
	printf("%s  \n",out);
}


void updVelOld(Cuerpo cuerpos[], int arraySize, double deltaT){
	int i;
	Cuerpo c;
	for (i=0;i<arraySize;i++){
		c = cuerpos[i];
		cuerpos[i].vel.x=cuerpos[i].vel.x+cuerpos[i].acel.x*deltaT;
		cuerpos[i].vel.y=cuerpos[i].vel.y+cuerpos[i].acel.y*deltaT;
		cuerpos[i].vel.z=cuerpos[i].vel.z+cuerpos[i].acel.z*deltaT;
	}
}
