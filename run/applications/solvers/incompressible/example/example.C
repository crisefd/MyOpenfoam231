/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"

    ///////Tensor operations for GeometricFields///////
    
    // U plus U
    //volVectorField UplusU("UplusU", U + U);
    
    //U minius U
    
    //volVectorField UminusU("UminusU", U - U);
    
    // U by p
    //volVectorField UbyP("UbyP", U * p);
    
    //U divide p
    //volVectorField UdividedP("UdividedP", U / p);
    
        // Magnitude of velocity field U
       //volScalarField magU("magU", mag(U));

       //Magnitude square of velocity field U
       //volScalarField magSqrU("magSqr", magSqr(U));
       
       // p power 3
       //volScalarField powP("powP", pow(p, 3));

        //Square of U
        //volSymmTensorField sqrU("sqrU", sqr(U));

        // U outer product U
        volTensorField UbyU("U*U", U * U);

        // U inner product U
        //volScalarField UdotU("UdotU", U & U);
        
        // U double inner product U
        //volScalarField UdotdotU("UdotdotU", UbyU && UbyU);
        
        // U cross product U
        //volVectorField UcrossU("UcrossU", U ^ U);

        //tr() Trace only works for spherical tensors

        //Min of U
        //dimensionedVector minU("minU", min(U));

        //Max of U
        //dimensionedVector maxU("maxU", max(U));

         //Transpose of UbyU?
        //volTensorField TransposeU("TransposeU", UbyU.T());

        //hodge dual of U
        //volTensorField hodgeDualU("hodgeDualU", *U);

    ///////////////////////////////////////////////

    //////Accessing the data in GeometricFields///////////////////
       /*label i = 0;
    //====Accessing the GeometricBoundaryField object of volField<Type> U (velocity field)
       volVectorField::GeometricBoundaryField boundaryFieldU(U.boundaryField());
       // Accessing the ith element fvPatchField<Type> of the boundary field of U
       const fvPatchVectorField & fvp = boundaryFieldU[i];
       label j = 0;
       // Accessing the jth element Type of the patch field of U
       Vector<double> vec1 = fvp[j];
       //Accessing component x of vector
       scalar vec1X = vec1.x();
       //creating a dimensioned<Type>  object from vector component
       dimensionedScalar dimScalar(
                   "dimScalar",
                   dimensionSet(0, 1, -1, 0, 0),
                   vec1X
                   );
    //====Accessing Internal Fields of volField<Type> object=======
        //Accessing the Field<Type> of GeometricField<Type> U
        Field<vector>internalFieldU(U.internalField());
        label k = 0;
        //Accessing the kth element Vector<Type> object of the internal Field
        Vector<double> vec2 = internalFieldU[k];
        //Accessing y component of vector
        scalar vec2Y = vec2.y();
        //Accessing the DimensionedInternalField<Type> object of GeometricField<Type> U
        volVectorField::DimensionedInternalField dimInternalFieldU(U.dimensionedInternalField());
        // Accessing dimensionSet object of the DimensionedInternalField object
        dimensionSet dimensions(dimInternalFieldU.dimensions());
        //Accessing the Mass dimension
        int mass = dimensions.MASS;
        //cout << "Mass dimension of DimensionedInternalField "<< mas << endl;
        //========Accessing the fvMesh object
        label inletI = mesh.boundaryMesh().findPatchID("inlet");
        scalar magInlet_1 = 0.0;
        scalar magInlet_2 = 0.0;
        const fvPatchVectorField & fvp_1 = U.mesh().C().boundaryField()[inletI];
        const fvPatchVectorField & fvp_2 = mesh.C().boundaryField()[inletI];
        if(fvp_1.size()){
            magInlet_1 = mag(fvp_1[0]);
        }
        if(fvp_2.size()){
            magInlet_1 = mag(fvp_2[0]);
        }
        Foam::Info << "mangInlet_1:" << magInlet_1 << " magInlet_2:" << magInlet_2 << Foam::endl;
        tensor t1(1, 2, 3, 4, 5, 6, 7, 8, 9);
        vector v1(-5, 5, 6);
        Foam::Info<<"v1.x():" << v1.x() <<" v1.component(x):"<<v1.component(vector::X)<< Foam::endl;
        Foam::Info<<"t1.xx():" << t1.xx() <<" t1.component(x):"<<t1.component(tensor::XX)<< Foam::endl;
        dimensionedScalar dim_scalar(
                    "dim_scalar",
                    dimensionSet(0, 1, -1, 0, 0),
                    t1.xx()
                    );


*/
///Finite volume calculus
   //gradient of scalar field p
   //vectorField gradP = fvc::grad(p);
    //volVectorField gradPhi = fvc::grad(phi);

   //sn gradient of scalar field p
    //surfaceScalarField snGradP(fvc::snGrad(p));

   //sn gradient of vector field U
   //surfaceVectorField snGradU(fvc::snGrad(U));

   //sn gradient of tensor field UbyU
    //surfaceTensorField snGradUbyU(fvc::snGrad(UbyU));

    //curl of U
    //volVectorField curlU(fvc::curl(U));

    //curl of p
    //volScalarField curlP(fvc::curl(p));



    //laplacian of scalar field p
    //volScalarField lapP(fvc ::laplacian(p));
     //volScalarField lapP_(fvc::laplacian(U, p));

     // Time deriative for scalar field p, vector field U and tensor field UbyU
      //volScalarField ddtP(fvc::ddt(p));
      //volScalarField ddtP_(fvc::ddt(p, p));
      //volVectorField ddtU(fvc::ddt(U));
      //volVectorField ddtU_(fvc::ddt(p, U));
      //volTensorField ddtUbyU(fvc::ddt(UbyU));
      //volTensorField ddtUbyU_(fvc::ddt(p, UbyU));

     //Second time derivative
      //volScalarField d2dt2P(fvc::d2dt2(p));
      //volScalarField d2dt2P_(fvc::d2dt2(p, p));
      //volVectorField d2dt2U(fvc::d2dt2(U));
      //volVectorField d2dt2U_(fvc::d2dt2(phi, U));
      //volTensorField d2dt2UbyU(fvc::d2dt2(UbyU));
      //volTensorField d2dt2UbyU_(fvc::d2dt2(p, UbyU));

     //Conective term(?)
      //volScalarField divP(fvc::div(phi, p));
      //volVectorField divU(fvc::div(phi, U));

      //volTensorField divUbyU(fvc::div(phi, UbyU));

     // divergent (?)
      //volScalarField divU_(fvc::div(U));
      //volVectorField divUbyU_(fvc::div(UbyU));
      //volScalarField divPhi(fvc::div(phi));


     //Source
      //volScalarField sourceP(fvc::Sp(2, p));
      // volScalarField sourceP_(fvc::Sp(p, p));
      //volVectorField sourceU(fvc::Sp(2, U));
      //volVectorField sourceU_(fvc::Sp(p, U));
       //volTensorField sourceUbyU(fvc::Sp(2, UbyU));
      //volTensorField sourceUbyU(fvc::Sp(p, UbyU));

      //volScalarField sourceP1(fvc::SuSp(p, p));
     //volVectorField sourceU1(fvc::SuSp(p, U));
      //volTensorField sourceUbyU1(fvc::SuSp(p, UbyU));

///Finite volume method

   //Time derivative of U
   //tmp<fvVectorMatrix> fvDdtU(fvm::ddt(U));
   //tmp<fvVectorMatrix> fvDdtU_(fvm::ddt(p, U));


   //Second time derivative of U
   //tmp<fvVectorMatrix> fvD2dt2U(fvm::d2dt2(U));
    //tmp<fvVectorMatrix> fvD2dt2U_(fvm::d2dt2(p, U));

   //Convective term(?) of phi and U
   //tmp<fvVectorMatrix> fvDivU(fvm::div(phi, U));

   //Laplacian
   //tmp<fvVectorMatrix> fvLaplacianU(fvm::laplacian(U));
   //Laplacian of p and U
   //tmp<fvVectorMatrix> fvLaplacianPU(fvm::laplacian(p, U));
   //Source of U
   //tmp<fvVectorMatrix> fvSpU(fvm::Sp(4, U));
   //tmp<fvVectorMatrix> fvSpU_(fvm::Sp(p, U));

   //Source depending of the sign
     //tmp<fvVectorMatrix> fvSuSpU_(fvm::SuSp(p, U));

///Operadores adicionales para fvMatrix

   //Sum with fvMatrix
     //tmp<fvVectorMatrix> fvSumExample_1(fvDdtU + U);
     //tmp<fvVectorMatrix> fvSumExample_2(fvDdtU + fvDdtU);

    //Substraction with fvMatrix
     //tmp<fvVectorMatrix> fvSubsExample_1(fvDdtU - U);
     //tmp<fvVectorMatrix> fvSubsExample_2(fvDdtU - fvDdtU);

   //Multiplication with fvMatrix
   //tmp<fvVectorMatrix> fvMultExample_1(p * fvDdtU );

///Methods for fvMatrix

   //Matrix scalar diagonal
   //scalarField diagonal = fvDdtU().D();
   //Matrix scalar upper side
   //scalarField upper = fvDdtU().upper();
   //Matrix scalar lower side
   //scalarField lower = fvDdtU().lower();
   //GeometricField central coefficients
   //volScalarField centralCoeff("AU",fvDdtU().A());
   //Matrix type diagonal
   //vectorField diagonal_ = fvDdtU().DD();
   //H operation source
   //volVectorField HoperationSource("HU",fvDdtU().H());
   //Matrix residual
   //vectorField residual = fvDdtU().residual();
   //face-flux field from the matrix
   //surfaceVectorField flux = fvDdtU().flux();

////////////////simple solver //////////////////////////////////////////////////////////
        simpleControl simple(mesh);
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
        Info<< "\nStarting time loop\n" << endl;
        while (simple.loop())
        {
            Info<< "Time = " << runTime.timeName() << nl << endl;
            // --- Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }
            turbulence->correct();
            runTime.write();
            Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                << nl << endl;
        }
        Info<< "End\n" << endl;

return 0;
}
