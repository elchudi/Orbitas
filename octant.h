#ifndef OCTANT
#define OCTANT

typedef struct tOctant{
	double xEnd;
	double yEnd;
	double zEnd;
	double xStart;
	double yStart;
	double zStart;
    struct tOctant *parent;
    struct tOctant *childs[8];
    Body bodies[1000];
    int isLeaf;
    int isRoot;
} Octant  ;

#endif  
