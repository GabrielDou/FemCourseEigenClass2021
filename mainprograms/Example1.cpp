#include <iostream>
#include "GeoMesh.h"
#include "GeoElementTemplate.h"
#include "Geom1d.h"
#include "Geom0d.h"
#include "VTKGeoMesh.h"
#include "CompMesh.h"
#include "Poisson.h"
#include "Analysis.h"
#include "PostProcessTemplate.h"

using namespace std;

int main(){
    
    // Creating a Geometric mesh
    int nElements = 64;
    int nNodes = nElements+1;

    GeoMesh *gmesh = new GeoMesh();   
    gmesh->SetDimension(1);
    gmesh->SetNumNodes(nNodes);
    gmesh->SetNumElements(nElements);

    VecDouble coord(3);
    coord[0] = 0.0;
    coord[1] = 0.0;
    coord[2] = 0.0;
    gmesh->Node(0).SetCo(coord);

    double deltaX = 1./nElements;
    double deltaY = 0.0;

    int matid = 1;

    for (int iEl = 0; iEl < nElements; iEl++)
    {
        coord[0] = deltaX*(iEl+1);
        coord[1] = deltaY*(iEl+1);
        gmesh->Node(iEl+1).SetCo(coord);

        VecInt nodes(2);
        nodes[0] = iEl;
        nodes[1] = iEl + 1;
        GeoElement *gel = new GeoElementTemplate<Geom1d>(nodes, matid, gmesh, iEl);
        gmesh->SetElement(iEl, gel);
    }

    //Inserindo os dois pontos da extremidade
    int matIdBC = 2;
    VecInt nIndex(1);
    //Primeiro ponto
    nIndex[0] = 0;
    GeoElement *point1 = new GeoElementTemplate<Geom0d>(nIndex, matIdBC, gmesh, nElements);
    gmesh->SetElement(nElements, point1);
    //Ultimo ponto
    nIndex[0] = nNodes-1;
    GeoElement *point2 = new GeoElementTemplate<Geom0d>(nIndex, matIdBC, gmesh, nElements+1);
    gmesh->SetElement(nElements+1, point2);

    gmesh->BuildConnectivity();

    // gmesh->Print(std::cout);

    VTKGeoMesh printer;
    printer.PrintGMeshVTK(gmesh,"geoMesh.vtk");
    
    //*************NEW*************
    //Cria uma malha computacional
    CompMesh cmesh(gmesh);

    //Estabelece um tipo de material
    MatrixDouble perm(3,3);
    perm.setIdentity();
    Poisson *mat = new Poisson(matid,perm);

    //Inserir o material de condição de contorno
    auto force = [](const VecDouble &x, VecDouble &res)
    {
        res[0] = x[0];
    };
    mat->SetForceFunction(force);
    MatrixDouble proj(1,1),val1(1,1),val2(1,1);
    proj.setZero();
    val1.setZero();
    val2.setZero();
    int bcType = 0;
    L2Projection *bc_point = new L2Projection(bcType,matIdBC,proj,val1,val2);
    std::vector<MathStatement *> mathvec = {0,mat,bc_point};
    cmesh.SetMathVec(mathvec);
    cmesh.SetDefaultOrder(2);
    
    //Insere o material na malha computacional e cria o espaço de aproximação
    cmesh.AutoBuild();

    //Análise
    Analysis an(&cmesh);
    an.RunSimulation();


    PostProcessTemplate<Poisson> postprocess;
    auto exact = [](const VecDouble &x, VecDouble &val, MatrixDouble &deriv)
    {
        val[0] = x[0]-sinh(x[0])/sinh(1.);
        deriv(0,0) = 1.-cosh(x[0])/sinh(1.);
    };
    postprocess.AppendVariable("Sol");
    postprocess.AppendVariable("DSol");

    postprocess.SetExact(exact);
    mat->SetExactSolution(exact);    

    // //Novo, mudar permeabilidade.
    VecDouble errvec;
    mat->SetExactSolution(exact);
    mat->SetDimension(1);
    an.PostProcessSolution("resultado64k2.vtk",postprocess);
    errvec = an.PostProcessError(std::cout, postprocess);


    return 0;
}


