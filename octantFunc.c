#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vector.h"
#include "body.h"
#include "octant.h"
#include "octantFunc.h"

/*Octant constructor !!! is this the way to do it? (i.e., if I need to create and array of poiters of structs in a function, should i do this or is ther such a thing as 'new struct' */
struct tOctant* getOctant (double xS, double yS, double zS, double xE, double yE, double zE, Octant *parent){
    Octant *toRet;
    toRet = malloc(sizeof(Octant));
    (*toRet).xStart = xS;
    (*toRet).yStart = yS;
    (*toRet).zStart = zS;
    (*toRet).xEnd = xE;
    (*toRet).yEnd = yE,
    (*toRet).zEnd = zE;
    (*toRet).parent = parent;
    return toRet;
}

Octant getRootOctant(Body bodies[], int qBodies){
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

    /*Correct way of linking arrays? !!!*/
    memcpy(&(c.bodies[0]),&(bodies[0]),qBodies);
    printf("qBodies : %d , body 0: %g , 1: %g\n",qBodies, c.bodies[0].masa, c.bodies[10].masa); 
    return c;
}

/*
void octantPartition(Body bodies[], int qBodies); 

void octantPartition(Body bodies[], int qBodies){
    Octant spaceT;
    spaceT = space(bodies, qBodies);

} 
*/

/*given an octant, it will populate the children field with 8 octants*/
void generateChildren(Octant *o){
    int i = 0, j = 0, k = 0, slices = 0, cont = 0;
    double xLen = 0, yLen = 0, zLen = 0;
    (*o).isLeaf = 0;
    slices = 2;
    if(slices != 0){
        xLen = ((*o).xEnd - (*o).xStart)/slices;
        yLen = ((*o).yEnd - (*o).yStart)/slices;
        zLen = ((*o).zEnd - (*o).zStart)/slices;
    }else{
        xLen = (*o).xEnd - (*o).xStart;
        yLen = (*o).yEnd - (*o).yStart;
        zLen = (*o).zEnd - (*o).zStart;
    }
    cont = 0;
    for(i = 0; i < slices; i++){
        for (j = 0; j < slices; j++){
            for (k = 0; k < slices; k++){
                
                /*
                aux.xStart = (*o).xStart + xLen*i;    
                aux.yStart = (*o).yStart + yLen*j;    
                aux.zStart = (*o).zStart + zLen*k;    
                aux.xEnd = (*o).xStart + xLen*(i+1);
                aux.yEnd = (*o).yStart + yLen*(j+1);
                aux.zEnd = (*o).zStart + zLen*(k+1);
                */
                (*o).childs[cont] = getOctant((*o).xStart + xLen*i, (*o).yStart + yLen*j, (*o).zStart + zLen*k, (*o).xStart + xLen*(i+1), (*o).yStart + yLen*(j+1), (*o).zStart + zLen*(k+1),o);
                cont++;
            }
        }
    } 

}

/*
void populateOctant(Octant *o){
    if((*o).isLeaf){
        
    
    }
    

}*/
