{
    volScalarField densityVariation
    (
        betaD*((rhoD - rhoC)/rho1)
    );

    surfaceScalarField R
    (
        fvc::interpolate(densityVariation)*phiVdj
    );


    scalar magDivRwA =
        mag
        (
            fvc::div(R)
        )().weightedAverage(rho*mesh.V()).value();

    scalar magDivRwADt =
        runTime.deltaTValue()*
        mag
        (
            fvc::div(R)
        )().weightedAverage(rho*mesh.V()).value();

    sumMagDivRwADt += magDivRwADt;


    scalar divRwA =
        (
            fvc::div(R)
        )().weightedAverage(rho*mesh.V()).value();

    scalar divRwADt =
        runTime.deltaTValue()*
        (
            fvc::div(R)
        )().weightedAverage(rho*mesh.V()).value();

    sumDivRwADt += divRwADt;


    Info<< "divR_calculation " << runTime.timeName()
        << " " << magDivRwA
        << " " << magDivRwADt
        << " " << sumMagDivRwADt
        << " " << divRwA
        << " " << divRwADt
        << " " << sumDivRwADt
        << endl;
}
