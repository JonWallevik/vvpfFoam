label nAlpha1Corr(readLabel(pimple.dict().lookup("nAlpha1Corr")));

label nAlpha1SubCycles(readLabel(pimple.dict().lookup("nAlpha1SubCycles")));


if (nAlpha1SubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum(0.0*rhoPhi);
    surfaceScalarField phiAlpha1Sum(0.0*phiAlpha1Global);

    for
    (
        subCycle<volScalarField> alpha1SubCycle(alpha1, nAlpha1SubCycles);
        !(++alpha1SubCycle).end();
    )
    {
        #ifdef ALPHA1RHO_SOLVE
        #include "alpha1EqnRho.H"
        #else
        #include "alpha1Eqn.H"
        #endif

        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;

        phiAlpha1Sum += (runTime.deltaT()/totalDeltaT)*phiAlpha1Global;
    }

    rhoPhi = rhoPhiSum;
    phiAlpha1Global = phiAlpha1Sum;
}
else
{
    #ifdef ALPHA1RHO_SOLVE
    #include "alpha1EqnRho.H"
    #else
    #include "alpha1Eqn.H"
    #endif
}

alpha1 == min(max(alpha1, scalar(0)), scalar(1));

rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
