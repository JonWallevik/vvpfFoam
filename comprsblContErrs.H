{
    // mag(div(phi)):
    scalar magWeightedAverageDivPhi =
        mag
        (
            fvc::div(phi)
        )().weightedAverage(rho*mesh.V()).value();

    // mag(dt*div(phi)):
    scalar magWeightedAverageDivPhiDt =
        runTime.deltaTValue()*
        mag
        (
            fvc::div(phi)
        )().weightedAverage(rho*mesh.V()).value();

    // sum(mag(div(phi)*dt)):
    sumMagWeightedAverageDivUdt += magWeightedAverageDivPhiDt;


    // div(phi):
    scalar weightedAverageDivPhi =
        (
            fvc::div(phi)
        )().weightedAverage(rho*mesh.V()).value();

    // dt*div(phi):
    scalar weightedAverageDivPhiDt =
        runTime.deltaTValue()*
        (
            fvc::div(phi)
        )().weightedAverage(rho*mesh.V()).value();

    // sum(div(phi)*dt):
    sumWeightedAverageDivUdt += weightedAverageDivPhiDt;


    Info<< "continuity_errors " << runTime.timeName()
        << " " << magWeightedAverageDivPhi
        << " " << magWeightedAverageDivPhiDt
        << " " << sumMagWeightedAverageDivUdt
        << " " << weightedAverageDivPhi
        << " " << weightedAverageDivPhiDt
        << " " << sumWeightedAverageDivUdt
        << endl;
}
