struct tOctant* getOctant (double xS, double yS, double zS, double xE, double yE, double zE, Octant *parent);
Octant getRootOctant(Body bodies[], int qBodies);
void generateChildren(Octant *o);
void populateOctant();
