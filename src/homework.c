#include "fem.h"

double geoSize(double x, double y) {
    femGeo *theGeometry = geoGetGeometry();
    double dc = theGeometry->dhCenter;
    double hc = theGeometry->hCenter;
    double forcePosx = theGeometry->forcePositionX;
    double forcePosy = theGeometry->forcePositionY;
    double forceR = theGeometry->forceRadius;
    double dt = theGeometry->dhTooth;
    double ht = theGeometry->hTooth;
    double h = theGeometry->h;
    double curvRatio = theGeometry->curvatureRatio;
    double curv = theGeometry->curvature;
    double rin = theGeometry->Rinner;
    double rout = theGeometry->Router;
    double toothLength = theGeometry->toothL;
    
    double distcenter = sqrt((y)*(y) + (x)*(x));
    double htemp = h;
    if(distcenter <= rin+dc)
    {
      double alpha = curvRatio * (2.0*M_PI/3.0);
      double beta = 2.0*asin(rin*curv*sin(alpha/2.0));
      double dist = rin*cos(alpha/2.0) + cos(beta/2.0)/curv;
      double distInner1 = sqrt((y-dist)*(y-dist) + (x)*(x));
      double distInner2 = sqrt((y-dist*cos(2.0*M_PI/3.0))*(y-dist*cos(2.0*M_PI/3.0)) + (x-dist*sin(2.0*M_PI/3.0))*(x-dist*sin(2.0*M_PI/3.0)));
      double distInner3 = sqrt((y-dist*cos(4.0*M_PI/3.0))*(y-dist*cos(4.0*M_PI/3.0)) + (x-dist*sin(4.0*M_PI/3.0))*(x-dist*sin(4.0*M_PI/3.0)));
      double distEdge = distcenter - rin;
  
      double tolerance = 1.0;
      if(distInner1 <= 1/curv+tolerance){
        double theta = asin(x/distInner1);
        if(theta < alpha/2.0 && theta > -alpha/2.0){
          distEdge = 1/curv - distInner1;
        }
      }
      if(distInner2 <= 1/curv+tolerance){
        double theta = M_PI/3.0 - atan((dist*sin(2.0*M_PI/3.0)-x)/(y-dist*cos(2.0*M_PI/3.0)));
        if(theta < alpha/2.0 && theta > -alpha/2.0){
          distEdge = 1/curv - distInner2;
        }
      }
      if(distInner3 <= 1/curv+tolerance){
        double theta = -M_PI/3.0 - atan((dist*sin(4.0*M_PI/3.0)-x)/(y-dist*cos(4.0*M_PI/3.0)));
        if(theta < alpha/2.0 && theta > -alpha/2.0){
          distEdge = 1/curv - distInner3;
        }
      }
      if(distEdge <= dc)
      {
        htemp = hc + 3*(h-hc)*((distEdge/dc)*(distEdge/dc)) + 2*(hc-h)*((distEdge/dc)*(distEdge/dc)*(distEdge/dc));
        h = fmin(htemp, h);
      }
    }
    
    htemp = h;
    double distForce = sqrt((y-forcePosy)*(y-forcePosy) + (x-forcePosx)*(x-forcePosx));
    if(distForce <= forceR)
    {
      double distEdge = rout+toothLength-distcenter;
      if(distEdge <= dt)
      {
        htemp = ht + 3*(h-hc)*((distEdge/dt)*(distEdge/dt)) + 2*(hc-h)*((distEdge/dt)*(distEdge/dt)*(distEdge/dt));
        h = fmin(htemp, h);
      }
  
    }
  
    return h;
  
  }

void geoMeshGenerate(double lc) {
    femGeo *theGeometry = geoGetGeometry();
    double ri = theGeometry->Rinner;
    double ro = theGeometry->Router;
    double curv = theGeometry->curvature;
    double curvRatio = theGeometry->curvatureRatio;
    int nTeeth = theGeometry->nTooths;
    double toothLength = theGeometry->toothL;
    double toothWidth = theGeometry->toothW;
  
    geoSetSizeCallback(geoSize);
  
    int ierr;
  
    // Create disk structure
    double alpha = curvRatio * (2.0*M_PI/3.0);
    if(ri*sin(alpha/2.0) > 1/curv)
    {
      printf("Error: curvRatio is too high\n");
      return;
    }
    double beta = 2.0*asin(ri*curv*sin(alpha/2.0));
    double dist = ri*cos(alpha/2.0) + cos(beta/2.0)/curv;
    int idOutDisk = gmshModelOccAddDisk(0,0,0,ro,ro,-1,NULL, (size_t)0 ,NULL,(size_t)0, &ierr);
    ErrorGmsh(ierr);
    int idInDisk = gmshModelOccAddDisk(0,0,0,ri,ri,-1,NULL,(size_t)0,NULL,(size_t)0, &ierr);
    ErrorGmsh(ierr);
    int idCurvDisk  = gmshModelOccAddDisk(0,dist,0,1/curv,1/curv,-1,NULL,(size_t)0,NULL,(size_t)0, &ierr);
    ErrorGmsh(ierr);
    int outDisk[] = {2, idOutDisk};
    int inDisk[] = {2, idInDisk};
    int curvDisk[] = {2, idCurvDisk};
  
    gmshModelOccCut(inDisk, 2,curvDisk, 2, NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr);
    ErrorGmsh(ierr);
    for(int k = 0; k < 2; k++)
    {
      gmshModelOccRotate(curvDisk, 2, 0, 0, 0, 0, 0, 1, 2*M_PI/3, &ierr);
      ErrorGmsh(ierr);
      gmshModelOccCut(inDisk, 2,curvDisk, 2, NULL,NULL,NULL,NULL,NULL,-1,1,0,&ierr);
      ErrorGmsh(ierr);
    }
    gmshModelOccRemove(curvDisk, 2, 0, &ierr);
    ErrorGmsh(ierr);
    int *holder;
    size_t holderdim;
    gmshModelOccCut(outDisk, 2,inDisk, 2, &holder,&holderdim,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
  
    // Create tooth structure
    int idTooth = gmshModelOccAddDisk(ro,0,0,toothLength,toothWidth,-1,NULL,(size_t)0,NULL,(size_t)0, &ierr);
    ErrorGmsh(ierr);
    int tooth[] = {2, idTooth};
    int idTorus = gmshModelOccAddDisk(0,0,0,2*ro,2*ro,-1,NULL,(size_t)0,NULL,(size_t)0, &ierr);
    ErrorGmsh(ierr);
    int torus[] = {2, idTorus};
    int idInTorus = gmshModelOccAddDisk(0,0,0,ro+(toothLength/1.25),ro+(toothLength/1.25),-1,NULL,(size_t)0,NULL,(size_t)0, &ierr);
    ErrorGmsh(ierr);
    int inTorus[] = {2, idInTorus};
    gmshModelOccCut(torus, 2,inTorus, 2, NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
    int *endTooth;
    size_t endToothdim;
    gmshModelOccCut(tooth, 2,torus, 2, &endTooth,&endToothdim,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
  
    for(int t = 0; t < nTeeth; t++)
    {
      // printf("t = %d\n", t);
      gmshModelOccRotate(endTooth, endToothdim, 0, 0, 0, 0, 0, 1, (2*M_PI/nTeeth), &ierr);
      ErrorGmsh(ierr);
      gmshModelOccSynchronize(&ierr);
      ErrorGmsh(ierr);
      int *temp;
      size_t tempdim;
      gmshModelOccFuse(holder, holderdim,endTooth, endToothdim, &temp,&tempdim,NULL,NULL,NULL,-1,1,0,&ierr);
      ErrorGmsh(ierr);
      holder = temp;
      holderdim = tempdim;
    }
    gmshModelOccRemove(endTooth, endToothdim, 0, &ierr);
    ErrorGmsh(ierr);
    // gmshModelOccRemove(holder, 2, 0, &ierr);
    // ErrorGmsh(ierr);
  
    gmshModelOccSynchronize(&ierr);
  
    // Use a frontal delaunay algorithm
    gmshOptionSetNumber("Mesh.Algorithm", 6, &ierr);
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);
  
    return;
}
  
void geoAssembleDomains(void){
  
    femGeo *theGeometry = geoGetGeometry();
    double rin = theGeometry->Rinner;
    double rout = theGeometry->Router;
    double limit = rin + (rout-rin)/2.0;
    double forceR = theGeometry->forceRadius;
    double forcePosx = theGeometry->forcePositionX;
    double forcePosy = theGeometry->forcePositionY;
    // printf("Limit: %f\n", limit);
  
    int domainAppartenance[theGeometry->nDomains];
    // printf("Number of domains: %d\n", theGeometry->nDomains);
    int innerNElem = 0;
    int freeNElem = 0;
    int forceNElem = 0;
  
    for(int i = 0; i < theGeometry->nDomains; i++)
    {
      domainAppartenance[i] = 0;
      femDomain *currentDomain = theGeometry->theDomains[i];
      int edge = currentDomain->elem[0];
      int nodeA = currentDomain->mesh->elem[2*edge];
      int nodeB = currentDomain->mesh->elem[2*edge+1];
      double xA = theGeometry->theNodes->X[nodeA];
      double yA = theGeometry->theNodes->Y[nodeA];
      double xB = theGeometry->theNodes->X[nodeB];
      double yB = theGeometry->theNodes->Y[nodeB];
      double x = xA + (xB-xA)/2.0;
      double y = yA + (yB-yA)/2.0;
      double dist = sqrt((x)*(x) + (y)*(y));
      double distForce = sqrt((x-forcePosx)*(x-forcePosx) + (y-forcePosy)*(y-forcePosy));
      // printf("Domain %i : %f\n", i, dist);
      if(dist < limit){
        domainAppartenance[i] = 1;
        innerNElem += currentDomain->nElem;
      }
      else if(distForce < forceR && xB < xA && fabs(yB-yA) < 0.35){
        domainAppartenance[i] = 2;
        forceNElem += currentDomain->nElem;
      }
      else{
        freeNElem += currentDomain->nElem;
      }
    }
    //newdomains
    femDomain *freeDomain = malloc(sizeof(femDomain));
    // freeDomain->name = "Free";
    freeDomain->mesh = theGeometry->theEdges;
    freeDomain->nElem = freeNElem;
    freeDomain->elem = malloc(sizeof(int) * freeNElem);
    femDomain *innerDomain = malloc(sizeof(femDomain));
    // freeDomain->name = "Inner";
    innerDomain->mesh = theGeometry->theEdges;
    innerDomain->nElem = innerNElem;
    innerDomain->elem = malloc(sizeof(int) * innerNElem);
    femDomain *forceDomain = malloc(sizeof(femDomain));
    forceDomain->mesh = theGeometry->theEdges;
    forceDomain->nElem = forceNElem;
    forceDomain->elem = malloc(sizeof(int) * forceNElem);
    int freeIndex = 0;
    int innerIndex = 0;
    int forceIndex = 0;
    for(int i = 0; i < theGeometry->nDomains; i++)
    {
      // printf("Domain %d: %d\n", i, domainAppartenance[i]);  
      femDomain *currentDomain = theGeometry->theDomains[i];
      if(domainAppartenance[i] == 1){
        for(int j = 0; j < currentDomain->nElem; j++){
          innerDomain->elem[innerIndex] = currentDomain->elem[j];
          innerIndex++;
        }
      }
      else if(domainAppartenance[i] == 2){
        for(int j = 0; j < currentDomain->nElem; j++){
          forceDomain->elem[forceIndex] = currentDomain->elem[j];
          forceIndex++;
        }
      }
      else{
        for(int j = 0; j < currentDomain->nElem; j++){
          freeDomain->elem[freeIndex] = currentDomain->elem[j];
          freeIndex++;
        }
      }
      currentDomain->nElem = 0;
      free(currentDomain->elem);
      currentDomain->mesh = NULL;
      free(currentDomain);
    }
    theGeometry->theDomains = realloc(theGeometry->theDomains, sizeof(femDomain *) * 3);
    theGeometry->theDomains[0] = freeDomain;
    theGeometry->theDomains[1] = innerDomain;
    theGeometry->theDomains[2] = forceDomain;
    theGeometry->nDomains = 3;
  
    return;
}

void femElasticityAssembleElements(femProblem *theProblem){
    // printf("Assembling elements...\n");
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    int nLocal = theMesh->nLocalNode;
    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    // printf("Starting assembly of elements...\n");
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                            dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                            dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                    A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                            dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                    A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                            dphidx[i] * c * dphidx[j]) * jac * weight; }}
             for (i = 0; i < theSpace->n; i++) {
                B[mapY[i]] -= phi[i] * g * rho * jac * weight; }}} 
}


void femElasticityAssembleNeumann(femProblem *theProblem){
    // printf("Assembling Neumann boundary conditions...\n");
    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->ruleEdge;
    femDiscrete    *theSpace = theProblem->spaceEdge;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    double x[2],y[2],phi[2];
    int iBnd,iElem,iInteg,iEdge,i,j,d,map[2],mapU[2];
    int nLocal = 2;
    double *B  = theSystem->B;

    for(iBnd=0; iBnd < theProblem->nBoundaryConditions; iBnd++){
        femBoundaryCondition *theCondition = theProblem->conditions[iBnd];
        femBoundaryType type = theCondition->type;
        double value = theCondition->value;

        if(type == DIRICHLET_X || type==DIRICHLET_Y){continue;}

        int shift = (type==NEUMANN_X) ? 0 : 1;
        // printf("Assembling Neumann boundary condition %d with %i elements\n", iBnd, theCondition->domain->nElem);

        //iterate over all elements
        for(iEdge = 0; iEdge < theCondition->domain->nElem; iEdge++){
            // printf("Assembling Neumann boundary condition %d on element %d\n", iBnd, iEdge);
            iElem = theCondition->domain->elem[iEdge];

            for(j=0; j < nLocal; j++) {
                map[j]  = theEdges->elem[iElem*nLocal+j];
                mapU[j] = 2*map[j] + shift;
                x[j]    = theNodes->X[map[j]];
                y[j]    = theNodes->Y[map[j]];
            }
                
            double jac = sqrt(pow(x[1]-x[0],2) + pow(y[1]-y[0],2))/2;

            for(iInteg=0; iInteg < theRule->n; iInteg++) {    
                double xsi    = theRule->xsi[iInteg];
                double weight = theRule->weight[iInteg];

                femDiscretePhi(theSpace,xsi,phi);

                for(i=0; i < theSpace->n; i++) {
                    B[mapU[i]] += phi[i] * value * jac * weight;
                }

        }

    }
}}

double **A_copy = NULL;
double *B_copy = NULL;

double *femElasticitySolve(femProblem *theProblem){

    femFullSystem  *theSystem = theProblem->system;
 
    femFullSystemInit(theSystem);
    // printf("After initialization\n");

    femElasticityAssembleElements(theProblem);
    // printf("After assembling elements\n");
    
    femElasticityAssembleNeumann(theProblem);
    // printf("After assembling Neumann conditions\n");

    int size = theSystem->size;

    if (A_copy == NULL)
    {
        A_copy = (double **) malloc(sizeof(double *) * size);
        for (int i = 0; i < size; i++) { A_copy[i] = (double *) malloc(sizeof(double) * size); }
    }
    if (B_copy == NULL) { B_copy = (double *) malloc(sizeof(double) * size); }

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            A_copy[i][j] = theSystem->A[i][j];
        }
        B_copy[i] = theSystem->B[i];
    }

    int *theConstrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < size; i++)
    {
        if (theConstrainedNodes[i] != -1)
        {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    femFullSystemEliminate(theSystem);
    memcpy(theProblem->soluce, theSystem->B, theSystem->size * sizeof(double));
    return theProblem->soluce;
}

double *femElasticityForces(femProblem *theProblem)
{
    double *residuals = theProblem->residuals;
    double *soluce    = theProblem->soluce;
    int size = theProblem->system->size;

    // Allocate memory for residuals if not already done
    if (residuals == NULL) { residuals = (double *) malloc(sizeof(double) * size); }

    // Initialize residuals to zero
    for (int i = 0; i < size; i++) { residuals[i] = 0.0; }

    /*
    Compute residuals: R = A * U - B where A and B are the system matrix
    and load vector before applying Dirichlet boundary conditions.
    */
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++) { residuals[i] += A_copy[i][j] * soluce[j]; }
        residuals[i] -= B_copy[i];
    }

    // Free memory allocated for the copy of the stiffness matrix A and the load vector B
    for (int i = 0; i < size; i++) { free(A_copy[i]); A_copy[i] = NULL;}
    free(A_copy); free(B_copy);
    A_copy = NULL; B_copy = NULL;

    // Return the residuals corresponding to the forces
    return residuals;
}

void renumberMesh(femGeo *theGeometry) {
    if (!theGeometry || !theGeometry->theNodes) {
        Error("Geometry or nodes not initialized");
    }

    int nNodes = theGeometry->theNodes->nNodes;

    //Build adjacency list
    femMesh *mesh = theGeometry->theElements ? theGeometry->theElements : theGeometry->theEdges;
    if (!mesh) {
        Error("No connectivity (elements/edges) available for RCM renumbering");
    }

    int nElem = mesh->nElem;
    int nLocal = mesh->nLocalNode;
    int *elem = mesh->elem;

    typedef struct {
        int count;
        int capacity;
        int *neighbors;
    } NodeAdjacency;

    NodeAdjacency *adjList = malloc(nNodes * sizeof(NodeAdjacency));
    if (!adjList) {
        Error("Memory allocation failed for adjacency list");
    }

    for (int i = 0; i < nNodes; i++) {
        adjList[i].count = 0;
        adjList[i].capacity = 4;
        adjList[i].neighbors = malloc(adjList[i].capacity * sizeof(int));
        if (!adjList[i].neighbors) {
            Error("Memory allocation failed for adjacency neighbors");
        }
    }

    for (int e = 0; e < nElem; e++) {
        for (int i = 0; i < nLocal; i++) {
            for (int j = i + 1; j < nLocal; j++) {
                int ni = elem[e * nLocal + i];
                int nj = elem[e * nLocal + j];

                // Add nj as a neighbor of ni
                int found = 0;
                for (int k = 0; k < adjList[ni].count; k++) {
                    if (adjList[ni].neighbors[k] == nj) {
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    if (adjList[ni].count == adjList[ni].capacity) {
                        adjList[ni].capacity *= 2;
                        adjList[ni].neighbors = realloc(adjList[ni].neighbors, adjList[ni].capacity * sizeof(int));
                    }
                    adjList[ni].neighbors[adjList[ni].count++] = nj;
                }

                // Add ni as a neighbor of nj
                found = 0;
                for (int k = 0; k < adjList[nj].count; k++) {
                    if (adjList[nj].neighbors[k] == ni) {
                        found = 1;
                        break;
                    }
                }
                if (!found) {
                    if (adjList[nj].count == adjList[nj].capacity) {
                        adjList[nj].capacity *= 2;
                        adjList[nj].neighbors = realloc(adjList[nj].neighbors, adjList[nj].capacity * sizeof(int));
                    }
                    adjList[nj].neighbors[adjList[nj].count++] = ni;
                }
            }
        }
    }

    //Compute RCM order
    int *rcmOrder = malloc(nNodes * sizeof(int));
    int *visited = calloc(nNodes, sizeof(int));
    if (!rcmOrder || !visited) {
        Error("Memory allocation failed for RCM order or visited array");
    }

    int rcmIndex = 0;
    for (int start = 0; start < nNodes; start++) {
        if (!visited[start]) {
            int *queue = malloc(nNodes * sizeof(int));
            if (!queue) {
                Error("Memory allocation failed for BFS queue");
            }

            int qStart = 0, qEnd = 0;
            queue[qEnd++] = start;
            visited[start] = 1;

            while (qStart < qEnd) {
                int current = queue[qStart++];
                rcmOrder[rcmIndex++] = current;

                int *neighbors = malloc(adjList[current].count * sizeof(int));
                if (!neighbors) {
                    Error("Memory allocation failed for neighbors");
                }

                int neighborCount = 0;
                for (int i = 0; i < adjList[current].count; i++) {
                    int neighbor = adjList[current].neighbors[i];
                    if (!visited[neighbor]) {
                        neighbors[neighborCount++] = neighbor;
                        visited[neighbor] = 1;
                    }
                }

                // Sort neighbors by degree
                for (int i = 0; i < neighborCount - 1; i++) {
                    for (int j = i + 1; j < neighborCount; j++) {
                        if (adjList[neighbors[i]].count > adjList[neighbors[j]].count) {
                            int temp = neighbors[i];
                            neighbors[i] = neighbors[j];
                            neighbors[j] = temp;
                        }
                    }
                }

                for (int i = 0; i < neighborCount; i++) {
                    queue[qEnd++] = neighbors[i];
                }

                free(neighbors);
            }

            free(queue);
        }
    }

    free(visited);

    // Reverse RCM order
    for (int i = 0; i < nNodes / 2; i++) {
        int temp = rcmOrder[i];
        rcmOrder[i] = rcmOrder[nNodes - i - 1];
        rcmOrder[nNodes - i - 1] = temp;
    }

    // Step 3: Create mapping from old to new indices
    int *mapping = malloc(nNodes * sizeof(int));
    if (!mapping) {
        Error("Memory allocation failed for mapping");
    }

    for (int i = 0; i < nNodes; i++) {
        mapping[rcmOrder[i]] = i;
    }

    //Update node coordinates
    double *newX = malloc(nNodes * sizeof(double));
    double *newY = malloc(nNodes * sizeof(double));
    if (!newX || !newY) {
        Error("Memory allocation failed for new coordinates");
    }

    for (int i = 0; i < nNodes; i++) {
        newX[i] = theGeometry->theNodes->X[rcmOrder[i]];
        newY[i] = theGeometry->theNodes->Y[rcmOrder[i]];
    }

    free(theGeometry->theNodes->X);
    free(theGeometry->theNodes->Y);
    theGeometry->theNodes->X = newX;
    theGeometry->theNodes->Y = newY;

    //Update connectivity
    if (theGeometry->theEdges) {
        femMesh *edges = theGeometry->theEdges;
        for (int i = 0; i < edges->nElem; i++) {
            for (int j = 0; j < edges->nLocalNode; j++) {
                edges->elem[i * edges->nLocalNode + j] = mapping[edges->elem[i * edges->nLocalNode + j]];
            }
        }
    }

    if (theGeometry->theElements) {
        femMesh *elements = theGeometry->theElements;
        for (int i = 0; i < elements->nElem; i++) {
            for (int j = 0; j < elements->nLocalNode; j++) {
                elements->elem[i * elements->nLocalNode + j] = mapping[elements->elem[i * elements->nLocalNode + j]];
            }
        }
    }

    //Free temporary resources
    for (int i = 0; i < nNodes; i++) {
        free(adjList[i].neighbors);
    }
    free(adjList);
    free(rcmOrder);
    free(mapping);
}



double* femElasticitySolveBandRCMK(femProblem* theProblem) {
    // Retrieve the full system
    femFullSystem* theSystem = theProblem->system;
    int n = theSystem->size;
    double **A = theSystem->A;
    double *B = theSystem->B;
    const double tolerance = 1e-12;

    // Initialize the system and assemble elements and Neumann conditions
    femFullSystemInit(theSystem);
    femElasticityAssembleElements(theProblem);
    femElasticityAssembleNeumann(theProblem);

    // Allocate memory for copies of A and B if not already done
    if (A_copy == NULL) {
        A_copy = (double **) malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            A_copy[i] = (double *) malloc(n * sizeof(double));
        }
    }
    if (B_copy == NULL) {
        B_copy = (double *) malloc(n * sizeof(double));
    }

    // Create copies of A and B for later use
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A_copy[i][j] = A[i][j];
        }
        B_copy[i] = B[i];
    }

    // Apply Dirichlet boundary conditions
    int *constrainedNodes = theProblem->constrainedNodes;
    for (int i = 0; i < n; i++) {
        if (constrainedNodes[i] != -1) {
            double value = theProblem->conditions[constrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem, i, value);
        }
    }

    // Calculate the effective bandwidth of the matrix
    int bandwidth = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (fabs(A[i][j]) > tolerance && (j - i) > bandwidth) {
                bandwidth = j - i;
            }
        }
    }
    printf("Band solver: bandwidth = %d\n", bandwidth);

    // Perform Gaussian elimination restricted to the bandwidth
    for (int k = 0; k < n; k++) {
        if (fabs(A[k][k]) <= tolerance) {
            printf("Pivot index %d, value %e\n", k, A[k][k]);
            Error("Cannot eliminate with such a pivot");
        }
        int rowLimit = (k + bandwidth + 1 < n) ? (k + bandwidth + 1) : n;
        for (int i = k + 1; i < rowLimit; i++) {
            double factor = A[i][k] / A[k][k];
            int colLimit = (k + bandwidth + 1 < n) ? (k + bandwidth + 1) : n;
            for (int j = k + 1; j < colLimit; j++) {
                A[i][j] -= factor * A[k][j];
            }
            B[i] -= factor * B[k];
        }
    }

    // Perform back substitution within the bandwidth
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        int colLimit = (i + bandwidth + 1 < n) ? (i + bandwidth + 1) : n;
        for (int j = i + 1; j < colLimit; j++) {
            sum += A[i][j] * B[j];
        }
        B[i] = (B[i] - sum) / A[i][i];
    }


    for(int i=0; i<n; i++){
        theProblem->soluce[i] = B[i];
    }

    return B;
}


void writeStressCSV(femProblem *theProblem, double *stressElem, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if(fp == NULL) {
        Error("Impossible d'ouvrir le fichier CSV en écriture.");
    }
    
    int nElem = theProblem->geometry->theElements->nElem;
    // Écrire l'en-tête du fichier CSV
    fprintf(fp, "Element,Sigma_x,Sigma_y,Tau_xy\n");
    
    // Parcourir tous les éléments et écrire les contraintes associées
    for (int i = 0; i < nElem; i++) {
        // Si stressElem est organisé par élément avec 3 composantes par élément,
        // alors l'élément i aura :
        // stressElem[3*i] = sigma_x, stressElem[3*i+1] = sigma_y, stressElem[3*i+2] = tau_xy.
        fprintf(fp, "%d,%e,%e,%e\n", i, stressElem[3*i], stressElem[3*i+1], stressElem[3*i+2]);
    }
    
    fclose(fp);
}
 
void writeStrainCSV(femProblem *theProblem, double *strain, const char *filename) {
    FILE *fp = fopen(filename, "w");
    if(fp == NULL) {
        Error("Impossible d'ouvrir le fichier CSV en écriture.");
    }
    
    int nElem = theProblem->geometry->theElements->nElem;
    // Écrire l'en-tête du fichier CSV
    fprintf(fp, "Element,epsilon_x,epsilon_y,epsilon_xy\n");
    
    // Parcourir tous les éléments et écrire les contraintes associées
    for (int i = 0; i < nElem; i++) {
        // Si stressElem est organisé par élément avec 3 composantes par élément,
        // alors l'élément i aura :
        // stressElem[3*i] = sigma_x, stressElem[3*i+1] = sigma_y, stressElem[3*i+2] = tau_xy.
        fprintf(fp, "%d,%e,%e,%e\n", i, strain[3*i], strain[3*i+1], strain[3*i+2]);
    }
    
    fclose(fp);
}


void calculateStrain(femProblem *theProblem, double **strainElem) {
    femMesh *mesh = theProblem->geometry->theElements;
    femNodes *nodes = theProblem->geometry->theNodes;
    femDiscrete *space = theProblem->space;
    femIntegration *rule = theProblem->rule;
    int nElem = mesh->nElem;
    int nLocal = mesh->nLocalNode;

    // Iterate over all elements
    for (int elem = 0; elem < nElem; elem++) {
        int map[4];
        double x[4], y[4];
        double u_local[8];

        // Map nodes and displacements for the current element
        for (int j = 0; j < nLocal; j++) {
            map[j] = mesh->elem[elem * nLocal + j];
            x[j] = nodes->X[map[j]];
            y[j] = nodes->Y[map[j]];
            u_local[2 * j] = theProblem->soluce[2 * map[j]];
            u_local[2 * j + 1] = theProblem->soluce[2 * map[j] + 1];
        }

        double strainSum[3] = {0.0, 0.0, 0.0};
        double totalWeight = 0.0;

        // Loop over integration points
        for (int ip = 0; ip < rule->n; ip++) {
            double xsi = rule->xsi[ip];
            double eta = (rule->eta != NULL) ? rule->eta[ip] : 0.0;
            double weight = rule->weight[ip];

            double phi[4], dphidxsi[4], dphideta[4];
            double dphidx[4], dphidy[4];

            // Compute shape functions and their derivatives
            femDiscretePhi2(space, xsi, eta, phi);
            femDiscreteDphi2(space, xsi, eta, dphidxsi, dphideta);

            // Compute Jacobian matrix
            double dxdxsi = 0.0, dxdeta = 0.0, dydxsi = 0.0, dydeta = 0.0;
            for (int i = 0; i < space->n; i++) {
                dxdxsi += x[i] * dphidxsi[i];
                dxdeta += x[i] * dphideta[i];
                dydxsi += y[i] * dphidxsi[i];
                dydeta += y[i] * dphideta[i];
            }
            double jacobian = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

            // Transform derivatives to global coordinates
            for (int i = 0; i < space->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jacobian;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jacobian;
            }

            // Compute displacement gradients
            double duxdx = 0.0, duydy = 0.0, duxdy = 0.0, duydx = 0.0;
            for (int i = 0; i < space->n; i++) {
                duxdx += dphidx[i] * u_local[2 * i];
                duydy += dphidy[i] * u_local[2 * i + 1];
                duxdy += dphidy[i] * u_local[2 * i];
                duydx += dphidx[i] * u_local[2 * i + 1];
            }

            // Compute strain components
            double eps_xx = duxdx;
            double eps_yy = duydy;
            double gamma_xy = (duxdy + duydx) / 2.0;

            // Accumulate weighted strain values
            strainSum[0] += eps_xx * jacobian * weight;
            strainSum[1] += eps_yy * jacobian * weight;
            strainSum[2] += gamma_xy * jacobian * weight;
            totalWeight += jacobian * weight;
        }

        // Store averaged strain values for the element
        (*strainElem)[elem * 3] = strainSum[0] / totalWeight;
        (*strainElem)[elem * 3 + 1] = strainSum[1] / totalWeight;
        (*strainElem)[elem * 3 + 2] = strainSum[2] / totalWeight;
    }
}
 
void calculateStress(femProblem *myProblem, double **stressElem, double *strainElem) {
    double E = myProblem->E;
    double nu = myProblem->nu;
    int nElem = myProblem->geometry->theElements->nElem;
    

    for (int iElem = 0; iElem < nElem; iElem++) {
        double eps_x = strainElem[3 * iElem];
        double eps_y = strainElem[3 * iElem + 1];
        double gamma_xy = strainElem[3 * iElem + 2];
        double sigma_x, sigma_y, tau_xy;
        if (myProblem->planarStrainStress == PLANAR_STRAIN) {
            sigma_x = (E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))) * eps_x +
                        (E * nu / ((1 + nu) * (1 - 2 * nu))) * eps_y;
            sigma_y = (E * nu / ((1 + nu) * (1 - 2 * nu))) * eps_x +
                        (E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))) * eps_y;
            tau_xy  = E / (2 * (1 + nu)) * gamma_xy;
        } else { // PLANAR_STRESS
            sigma_x = E / (1 - nu * nu) * (eps_x + nu * eps_y);
            sigma_y = E / (1 - nu * nu) * (eps_y + nu * eps_x);
            tau_xy  = E / (2 * (1 + nu)) * gamma_xy;
        }
        (*stressElem)[3 * iElem]     = sigma_x;
        (*stressElem)[3 * iElem + 1] = sigma_y;
        (*stressElem)[3 * iElem + 2] = tau_xy;
    }
}
