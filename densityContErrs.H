// See also alphaEqnsSubCycle.H in interMixingFoam:
{
    volScalarField rhoNew(rho);

    solve(fvm::ddt(rho) + fvc::div(rhoPhi));

    scalar maxDensityError =
        max
        (
            (mag(rho - rhoNew)/mag(rhoNew))().internalField()
        );

    scalar averageDensityError =
        mag
        (
            (rho - rhoNew)/mag(rhoNew)
        )().average().value();

    Info<< "density_errors " << runTime.timeName() << " "
        << maxDensityError << " "
        << averageDensityError
        << endl;

    rho == rhoNew;
}
