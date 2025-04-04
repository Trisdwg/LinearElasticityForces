#include "fem.h"


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
    printf("Assembling Neumann boundary conditions...\n");
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
    printf("After initialization\n");

    femElasticityAssembleElements(theProblem);
    printf("After assembling elements\n");
    
    femElasticityAssembleNeumann(theProblem);
    printf("After assembling Neumann conditions\n");

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

    return B;
}
