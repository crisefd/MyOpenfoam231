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
/*
///////Tensor operations for GeometricFields///////
    // Magnitude of velocity field U
   volScalarField magU("magU", mag(U));
   //Magnitude square of velocity field U
   volScalarField magSqrU("magSqr", magSqr(U));
///////////////////////////////////////////////

//////Accessing the data in GeometricFields///////////////////
   label i = 0;
//****Accessing the GeometricBoundaryField object of volField<Type> U (velocity field)
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

//****Accessing Internal Fields of volField<Type> object*******
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
*/
    //****Accessing the fvMesh object
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
   return 0;
}
