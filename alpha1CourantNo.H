scalar maxAlpha1Co
(
    readScalar(runTime.controlDict().lookup("maxAlpha1Co"))
);

scalar alpha1CoNum = 0.0;
scalar meanAlpha1CoNum = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        pos(alpha1 - 0.01)*pos(0.99 - alpha1)
       *fvc::surfaceSum(mag(phi))().internalField()
    );

    alpha1CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanAlpha1CoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "alpha1 Courant number mean: " << meanAlpha1CoNum
    << " max: " << alpha1CoNum << endl;
