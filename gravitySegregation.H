#include "correctViscosity.H"

tmp<volVectorField> VsGR =
    mag(alpha1)*(pow(Da, 2.0)*g*(rhoD - rhoC))/(18.0*(etaEff/1.45));


#ifdef GRAVITY_SEGREGATION
VdjGR = 
    slowDown2
    (
        alphaD, // alphaD, or betaD, depending on user preference!
        alphaDMIN,
        alphaDMAX
    )*(1.0*VsGR);
#else
VdjGR = zeroVelocity;
#endif


forAll(alpha1.internalField(), celli)
{
    if
    (
        alpha1[celli] > lowerCrit.value()
     && alpha1[celli] < upperCrit.value()
     && alphaD[celli] > criteriaD.value()
    )
    {
        VdjGR[celli] = (1.0 - alpha1[celli])*0.2100*interfaceNormal[celli];
    }
    else if (alpha1[celli] <= lowerCrit.value())
    {
        VdjGR[celli] = (1.0 - alpha1[celli])*0.0306*g.value();
    }
}


forAll(mesh.boundary(), patchi)
{
    VdjGR.boundaryField()[patchi] == vector::zero;

    forAll(alpha1.boundaryField()[patchi], facei)
    {
        if (alpha1.boundaryField()[patchi][facei] < lowerCrit.value())
        {
            VdjGR.boundaryField()[patchi][facei] =
                (1.0 - alpha1.boundaryField()[patchi][facei])*0.0306*g.value();
        }
    }
}


VdjGR.correctBoundaryConditions();

phiVdjGR = fvc::interpolate(VdjGR) & mesh.Sf();
