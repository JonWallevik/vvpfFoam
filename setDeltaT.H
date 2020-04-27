if (adjustTimeStep)
{
    scalar maxDeltaTFact1 =
        min(maxCo/(CoNum + SMALL), maxAlpha1Co/(alpha1CoNum + SMALL));

    scalar maxDeltaTFactD =
        min(maxCo/(CoNum + SMALL), maxAlphaDCo/(alphaDCoNum + SMALL));

    scalar maxDeltaTFact = min(maxDeltaTFact1, maxDeltaTFactD);

    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact*runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
