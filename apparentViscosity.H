/*---------------------------------------------------------------------------*\
Example of a linear weight function that can be used to increase or
decrease the apparent viscosity eta_1. For example,...
     If varPhi = varPhiMAX = e.g. 0.4 => WF1 = WF1_MAX = e.g. 1.6;
     If varPhi = varPhiMIN = e.g. 0.0 => WF1 = WF1_MIN = e.g. 0.4;
     If varPhi =             e.g. 0.2 => WF1 =           e.g. 1.0;
\*---------------------------------------------------------------------------*/

volScalarField weightFunction1
(
    const volScalarField& varPhi,
    const dimensionedScalar& varPhiMIN,
    const dimensionedScalar& varPhiMAX
)
{
    dimensionedScalar WF1_MAX
    (
        "WF1_MAX",
        dimensionSet(0,0,0,0,0,0,0),
        scalar(1.6)
    );

    dimensionedScalar WF1_MIN
    (
        "WF1_MIN",
        dimensionSet(0,0,0,0,0,0,0),
        scalar(0.4)
    );

    dimensionedScalar slope("slope", dimensionSet(0,0,0,0,0,0,0), scalar(1));
    slope = (WF1_MAX - WF1_MIN)/(varPhiMAX - varPhiMIN);

    tmp<volScalarField> weight
    (
        slope*mag(varPhi - varPhiMIN) + WF1_MIN
    );

    return weight();
}

/*---------------------------------------------------------------------------*\
Different types of apparent viscosity functions. The appropriate one is chosen
by uncommenting the corresponding function name in the return segment.
                                    - * -
The term "varPhi" (c.f. $\varphi$ in LaTex) designates volume fraction of the
dispersed phase. The term varPhi can either be alphaD or betaD, depending on
user preference. Modification is made in the source code correctViscosity.H.
See also particleMigration.H and gravitySegregation.H.
\*---------------------------------------------------------------------------*/

volScalarField apparentViscosity
(
    const dimensionedScalar& mu,
    const dimensionedScalar& tau0,
    const dimensionedScalar& delta,
    const dimensionedScalar& eta2,
    const volScalarField& shearRate,
    const volScalarField& varPhi,
    const dimensionedScalar& varPhiMIN,
    const dimensionedScalar& varPhiMAX,
    const volScalarField& alpha1
)
{
    // A linear weight function:
    volScalarField WF1 =
        weightFunction1
        (
            varPhi,
            varPhiMIN,
            varPhiMAX
        );


    // Plastic viscosity mu:
    const dimensionedScalar mu_Bi("mu_Bi", dimPressure*dimTime, scalar(50.0));
    // Yield stress tau0:
    const dimensionedScalar tau0_Bi("tau0_Bi", dimPressure, scalar(10.0));

    // Bingham model:
    tmp<volScalarField> viscous_1
    (
        mag(alpha1)*
        (
            mu_Bi + tau0_Bi/(shearRate + delta)
        )
      + mag(scalar(1) - alpha1)*eta2
    );

    // Bingham model mixed with weightFunction1:
    tmp<volScalarField> viscous_2
    (
        WF1*mag(alpha1)*
        (
            mu_Bi + tau0_Bi/(shearRate + delta)
        )
      + mag(scalar(1) - alpha1)*eta2
    );

    // Herschel-Bulkley model mixed with weightFunction1:
    const dimensionedScalar tOne
    (
        "tOne",
        dimensionSet(0,0,1,0,0,0,0),
        scalar(1.0)
    );
    // Consistency index n:
    const scalar n_HB = 1.14;
    // Consistency factor K:
    const dimensionedScalar K_HB("K_HB", dimPressure*dimTime, scalar(40.7));
    // Yield stress tau0:
    const dimensionedScalar tau0_HB("tau0_HB", dimPressure, scalar(16.5));

    tmp<volScalarField> viscous_3
    (
        WF1*mag(alpha1)*
        (
            K_HB*pow(shearRate*tOne,n_HB-1.0) + tau0_HB/(shearRate + delta)
        )
      + mag(scalar(1) - alpha1)*eta2
    );


    // Modified Bingham model mixed with weightFunction1:
    // Plastic viscosity mu:
    const dimensionedScalar mu_mBi("mu_mBi", dimPressure*dimTime, scalar(82.6));
    // Second order term c:
    const dimensionedScalar c_mBi
    (
        "c_mBi",
        dimPressure*dimTime*dimTime,
        scalar(1.5)
    );
    // Yield stress tau0:
    const dimensionedScalar tau0_mBi("tau0_mBi", dimPressure, scalar(23.7));

    // Apparent viscosity with modified Bingham fluid and weightFunction1:
    tmp<volScalarField> viscous_4
    (
        WF1*mag(alpha1)*
        (
            mu_mBi + c_mBi*shearRate + tau0_mBi/(shearRate + delta)
        )
      + mag(scalar(1) - alpha1)*eta2
    );


    // Maximum packing fraction:
    const dimensionedScalar varPhiM("varPhiM", dimless, scalar(0.55));

    // Intrinsic viscosity:
    const dimensionedScalar etaInVi("etaInVi", dimless, scalar(3.25));

    // Model of (A) Krieger and Dougherty (B) Chateau, Ovarlez and Trung:
    tmp<volScalarField> viscous_5
    (
        mag(alpha1)*
        (
            mu*pow(mag(scalar(1) - varPhi/varPhiM), -etaInVi*varPhiM)
          + tau0*sqrt
            (
                mag(scalar(1) - varPhi)*pow(mag(scalar(1) - varPhi/varPhiM),
               -2.5*varPhiM)
            // -etaInVi*varPhiM)
            )/(shearRate + delta)
        )
      + mag(scalar(1) - alpha1)*eta2
    );


    // Model of Krieger and Dougherty:
    tmp<volScalarField> viscous_6
    (
        mag(alpha1)*
        (
            // mu*pow(mag(scalar(1) - varPhi/varPhiM), -etaInVi*varPhiM)
            mu*pow(mag(scalar(1) - varPhi/scalar(0.68)),scalar(-1.82))
        )
      + mag(scalar(1) - alpha1)*eta2
    );


    return viscous_5();
    // return viscous_6();
}
