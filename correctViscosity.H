// The apparent viscosity is implemented in the source apparentViscosity.H
etaEff =
    apparentViscosity
    (
        mu,
        tau0,
        delta,
        eta2,
        shearRate,
        betaD,     // alphaD, or betaD, depending on user preference!
        alphaDMIN,
        alphaDMAX,
        alpha1
    );
