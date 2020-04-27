label nAlphaDCorr(readLabel(pimple.dict().lookup("nAlphaDCorr")));

label nAlphaDSubCycles(readLabel(pimple.dict().lookup("nAlphaDSubCycles")));


if (nAlphaDSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiAlphaDSum(0.0*rhoPhiAlphaD);

    for
    (
        subCycle<volScalarField> alphaDSubCycle(alphaD, nAlphaDSubCycles);
        !(++alphaDSubCycle).end();
    )
    {
        #include "alphaDEqn.H"
        rhoPhiAlphaDSum += (runTime.deltaT()/totalDeltaT)*rhoPhiAlphaD;
    }

    rhoPhiAlphaD = rhoPhiAlphaDSum;
}
else
{
    #include "alphaDEqn.H"
}


alphaD == min(max(alphaD, alphaDMIN), alphaDMAX);


rho1 == (mag(alpha1 - alphaD)*rhoC + alphaD*rhoD)/(alpha1 + delta1);
rho1 == min(max(rho1, rho1MIN), rho1MAX);

rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;

surfaceScalarField rho1f
(
    "rho1f",
    fvc::interpolate(rho1)
);

rhoPhi = phiAlpha1Global*(rho1f - rho2) + phi*rho2;
