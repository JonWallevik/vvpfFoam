#include "correctViscosity.H"


volVectorField gradShearRate(fvc::grad(shearRateAlpha1));
// volVectorField gradShearRate(fvc::grad(shearRate));

// ---
// const dimensionedScalar Ksr("Ksr", dimLength*dimLength, scalar(0.8e-5));
// tmp<volVectorField> VsSR = -mag(alpha1)*Ksr*gradShearRate;
// ---

const dimensionedScalar Kc("Kc", dimless, scalar(0.41));
const dimensionedScalar Keta("Keta", dimless, scalar(0.62));
const dimensionedScalar a = Da/2.0;

tmp<volVectorField> VsSR =
    -mag(alpha1)*
    (
        //   Kc*pow(a, 2.0)*fvc::grad(shearRate*alphaD)
        // + Keta*pow(a, 2.0)*shearRate*alphaD*(gradEtaEff/(etaEff + deltaEta))
        Kc*pow(a, 2.0)*fvc::grad(shearRateAlpha1*alphaD)
      + Keta*pow(a, 2.0)*shearRateAlpha1*alphaD
      * (gradEtaEff/(etaEff + deltaEta))
    );


#ifdef PARTICLE_MIGRATION
VdjSR =
    slowDown4
    (
        alphaD, // alphaD, or betaD, depending on user preference!
        alphaDMIN,
        alphaDMAX
    )*(1.0*VsSR);
#else
VdjSR = zeroVelocity;
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
        VdjSR[celli] = (1.0 - alpha1[celli])*0.2100*interfaceNormal[celli];
    }
    else if (alpha1[celli] <= lowerCrit.value())
    {
        VdjSR[celli] = (1.0 - alpha1[celli])*0.0306*g.value();
    }
}


forAll(mesh.boundary(), patchi)
{
    VdjSR.boundaryField()[patchi] == vector::zero;

    /*
    forAll(alpha1.boundaryField()[patchi], facei)
    {
        if (alpha1.boundaryField()[patchi][facei] < lowerCrit.value())
        {
            VdjSR.boundaryField()[patchi][facei] =
                (1.0 - alpha1.boundaryField()[patchi][facei])*0.0306*g.value();
        }
    }
    */
}


VdjSR.correctBoundaryConditions();

phiVdjSR = fvc::interpolate(VdjSR) & mesh.Sf();
