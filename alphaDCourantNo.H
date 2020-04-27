scalar maxAlphaDCo
(
    readScalar(runTime.controlDict().lookup("maxAlphaDCo"))
);

// volScalarField rho1RatioVol(rhoC/rho1);
// surfaceScalarField rho1Ratio(fvc::interpolate(rho1RatioVol));

scalar alphaDCoNum = 0.0;
scalar meanAlphaDCoNum = 0.0;

if (mesh.nInternalFaces())
{
    scalarField sumPhi
    (
        pos(alphaD - 0.01)*pos(0.99 - alphaD)
       *fvc::surfaceSum(mag(phi))().internalField()
    // *fvc::surfaceSum(mag(phi + rho1Ratio*phiVdj))().internalField()
    );

    alphaDCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

    meanAlphaDCoNum =
        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "alphaD Courant number mean: " << meanAlphaDCoNum
    << " max: " << alphaDCoNum << endl;
