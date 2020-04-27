{
    word alphaDScheme("div(phi,alphaD)");
    word phiVdjScheme("div(phiVdj,alphaD)");

    word alphaDrScheme("div(phirb,alphaD)");
    word alpha1rScheme("div(phirb,alpha1)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(cAlpha*phic, max(phic));
    surfaceScalarField phir(phic*nHatf);


    volScalarField rho1RatioVol
    (
        "rho1RatioVol",
        rhoC/rho1
    );

    surfaceScalarField rho1Ratio
    (
        "rho1Ratio",
        fvc::interpolate(rho1RatioVol)
    );

    for (int aCorr=0; aCorr<nAlphaDCorr; aCorr++)
    {
        surfaceScalarField phiAlphaD
        (
            fvc::flux
            (
                phi,
                alphaD,
                alphaDScheme
            )
          + fvc::flux
            (
                rho1Ratio*phiVdj,
                alphaD,
                phiVdjScheme
            )
          + fvc::flux
            (
                fvc::flux
                (
                    -fvc::flux(-phir, scalar(1) - alpha1, alpha1rScheme),
                    alpha1,
                    alpha1rScheme
                ),
                scalar(1.00/alphaDMAX.value())*alphaD,alphaDrScheme
            )
        );


        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar
            (
                "Sp",
                dimensionSet(0,0,-1,0,0,0,0),
                0.0
            )
        );

        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            mesh,
            dimensionedScalar
            (
                "Su",
                dimensionSet(0,0,-1,0,0,0,0),
                0.0
            )
        );


        #if defined(ALPHAD_EXPLICIT_SOLVE)
        MULES::explicitSolve
        (
            geometricOneField(),
            alphaD,
            phi,
            phiAlphaD,
            Sp,
            Su,
            alphaDMAX.value(),
            alphaDMIN.value()
        );
        #elif defined(ALPHAD_IMPLICIT_SOLVE)
        MULES::implicitSolve
        (
            geometricOneField(),
            alphaD,
            phi,
            phiAlphaD,
            Sp,
            Su,
            alphaDMAX.value(),
            alphaDMIN.value()
        );
        #elif defined(ALPHAD_MATRIX_SOLVE)
        // Note, the following method is highly non-conservative and
        // not recommended! Use only for testing/experimentation!
        surfaceScalarField interfaceD
        (
            "interfaceD",
            phir*fvc::interpolate
            (
                alpha1*(scalar(1) - alpha1)*scalar(1.00/alphaDMAX.value())
            )
        );

        fvScalarMatrix alphaDsysEqn
        (
            fvm::ddt(alphaD)
          + fvm::div(phi, alphaD, alphaDScheme)
          + fvm::div(rho1Ratio*phiVdj, alphaD, phiVdjScheme)
          + fvc::div(interfaceD, alphaD, alphaDrScheme)
       // + fvc::div(interfaceD, alphaD, alpha1rScheme)
         ==
            Sp*alphaD + Su
        );

        alphaDsysEqn.relax();
        alphaDsysEqn.solve();
        #endif


        rho1 == (mag(alpha1 - alphaD)*rhoC + alphaD*rhoD)/(alpha1 + delta1);
        rho1 == min(max(rho1, rho1MIN), rho1MAX);

        rhoPhiAlphaD = phiAlphaD*fvc::interpolate(rho1);
    }

    Info<< "Dispersed phase (alphaD) volume fraction = "
        << alphaD.weightedAverage(mesh.V()).value()
        << "  Min(alphaD) = " << min(alphaD).value()
        << "  Max(alphaD) = " << max(alphaD).value()
        << endl;
}
