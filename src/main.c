/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *  Calcul des densités de force aux noeuds contraints
 *
 *  Copyright (C) 2024 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"

double fun(double x, double y) 
{
    return 1;
}

int main(int argc, char *argv[])
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    X : Horizontal residuals for unconstrained equations \n");
    printf("    Y : Horizontal residuals for unconstrained equations \n");
    printf("    N : Next domain highlighted\n\n\n");

    char *inputfile = "../data/mesh2,7k.txt";
    char *problemfile = "../data/defaultProblem.txt";
    int useGaussSolver = 0;
    femElasticCase cas = PLANAR_STRESS;

    for(int i =0; i < argc; ++i){
        if(strncasecmp(argv[i], "-in",3)==0)
            inputfile = argv[i+1];
        if(strncasecmp(argv[i], "-gauss",5)==0)
            useGaussSolver = 1;
        if(strncasecmp(argv[i], "-p", 2)==0)
            problemfile = argv[i+1];
        if(strncasecmp(argv[i], "-axisym", 7)==0)
            cas = AXISYM;
        if(strncasecmp(argv[i], "-defoplane", 13)==0)
            cas = PLANAR_STRAIN;
    }

    printf("Input file: %s\n", inputfile);
    printf("Solver: %s\n", useGaussSolver == 1 ? "Gauss" : "BandRCMK");
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    // theGeometry->elementType = FEM_TRIANGLE;
    // 
    // Paramètres de la géométrie
    // double ri = 18.85;
    // double ro = 25.0;
    // double curv = 1.0/(1.0*ri);
    // double curvRatio = 0.4583333333333333;
    // int nTeeth = 32;
    // double toothLength = 3.0;
    // double toothWidth = 1.5;
    // theGeometry->Rinner = ri;
    // theGeometry->Router = ro;
    // theGeometry->curvature = curv;
    // theGeometry->curvatureRatio = curvRatio;
    // theGeometry->toothL = toothLength;
    // theGeometry->toothW = toothWidth;
    // theGeometry->nTooths = nTeeth;
    // double dc = 3.5;
    // theGeometry->dhCenter = dc;
    // double hc = 0.5;
    // theGeometry->hCenter = hc;
    // double forcePosx = 50.0;
    // theGeometry->forcePositionX = forcePosx;
    // double forcePosy = 0.0;
    // theGeometry->forcePositionY = forcePosy;
    // double forceR = 30.0;
    // theGeometry->forceRadius = forceR;
    // double dt = 10.0;
    // theGeometry->dhTooth = dt;
    // double ht = 0.5;
    // theGeometry->hTooth = ht;
    // double h = 5.0;
    // theGeometry->h = h;

    // geoFinalize();
    // geoInitialize();
    geoMeshRead(inputfile);
    
        
//
//  -2- Creation probleme 
//

    double E, nu, rho, g;
    double forceIntensity;
    FILE *file = fopen(problemfile, "r");
    ErrorScan(fscanf(file, "Young modulus: %le \n", &E));
    ErrorScan(fscanf(file, "Density: %le \n", &rho));
    ErrorScan(fscanf(file, "Poisson ratio: %le \n", &nu));
    ErrorScan(fscanf(file, "Force intensity: %le \n", &forceIntensity));
    ErrorScan(fscanf(file, "Gravity: %le\n", &g));
    fclose(file);
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,cas);
    renumberMesh(theProblem->geometry);
    femElasticityAddBoundaryCondition(theProblem, "Inner", DIRICHLET_X, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Inner", DIRICHLET_Y, 0.0);
    femElasticityAddBoundaryCondition(theProblem, "Force", NEUMANN_Y, forceIntensity);
    femElasticityPrint(theProblem);

//
//  -3- Resolution du probleme et calcul des forces
//

    double *theSoluce = (useGaussSolver == 1) ? femElasticitySolve(theProblem) : femElasticitySolveBandRCMK(theProblem);
    double *theForces = femElasticityForces(theProblem);
    double area = femElasticityIntegrate(theProblem, fun);   


    //ecrivons dans un csv

    double *stressElem;
    double *strainElem;
    if (theProblem->planarStrainStress == AXISYM) {
        stressElem = malloc(4 * theProblem->geometry->theElements->nElem * sizeof(double));
        strainElem = malloc(4 * theProblem->geometry->theElements->nElem * sizeof(double));
    } else {
        stressElem = malloc(3 * theProblem->geometry->theElements->nElem * sizeof(double));
        strainElem = malloc(3 * theProblem->geometry->theElements->nElem * sizeof(double));
    }

    calculateStrain(theProblem, &strainElem);
    calculateStress(theProblem, &stressElem, strainElem);

    writeStressCSV(theProblem, stressElem, "stress.csv");
    writeStrainCSV(theProblem, strainElem, "strain.csv");
    //free(stressElem);
   
//
//  -4- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//

    
    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1e5;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    double *forcesX = malloc(theNodes->nNodes * sizeof(double));
    double *forcesY = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]);
        forcesX[i] = theForces[2*i+0];
        forcesY[i] = theForces[2*i+1]; }
  
    double hMin = femMin(normDisplacement,theNodes->nNodes);  
    double hMax = femMax(normDisplacement,theNodes->nNodes);  
    printf(" ==== Minimum displacement          : %14.7e [m] \n",hMin);
    printf(" ==== Maximum displacement          : %14.7e [m] \n",hMax);

//
//  -5- Calcul de la force globaleresultante
//

    double theGlobalForce[2] = {0, 0};
    for (int i=0; i<theProblem->geometry->theNodes->nNodes; i++) {
        theGlobalForce[0] += theForces[2*i+0];
        theGlobalForce[1] += theForces[2*i+1]; }
    printf(" ==== Global horizontal force       : %14.7e [N] \n",theGlobalForce[0]);
    printf(" ==== Global vertical force         : %14.7e [N] \n",theGlobalForce[1]);
    printf(" ==== Weight                        : %14.7e [N] \n", area * rho * g);

//
//  -6- Visualisation du maillage
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Recovering forces on constrained nodes");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'X') == GLFW_PRESS) { mode = 2;}
        if (glfwGetKey(window,'Y') == GLFW_PRESS) { mode = 3;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        if (t-told > 0.5) {freezingButton = FALSE; }
        
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 2) {
            glfemPlotField(theGeometry->theElements,forcesX);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 3) {
            glfemPlotField(theGeometry->theElements,forcesY);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);
    free(forcesX);
    free(forcesY);
    femElasticityFree(theProblem) ; 
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 
