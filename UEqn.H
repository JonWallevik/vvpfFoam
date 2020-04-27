    volScalarField alphaDrho1Ratio
    (
        "alphaDrho1Ratio",
        (betaD/(scalar(1) - betaD))*((rhoC*rhoD)/rho1)
    );

    rhoPhiVdj = fvc::interpolate(alphaDrho1Ratio)*phiVdj;

    // For single reference frame calculations:
    #ifdef SINGLE_REFERENCE_FRAME
    volVectorField Fcoriolis = rho*(2*Omega ^ U);
    volVectorField Fcentrifugal = rho*(Omega ^ (Omega ^ mesh.C()));
    // volVectorField Feuler = rho*(dOmegaDt ^ mesh.C());
    #endif

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      + fvc::div(rhoPhiVdj, Vdj, "div(rhoPhiVdj,Vdj)")
      - fvm::laplacian(etaEff, U)
      - fvc::div(etaEff*dev(T(fvc::grad(U))))
        #ifdef SINGLE_REFERENCE_FRAME
      + Fcoriolis
      + Fcentrifugal
   // + Feuler
        #endif
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(sigmaK)*fvc::snGrad(alpha1)
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
