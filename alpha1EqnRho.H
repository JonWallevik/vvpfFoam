{
    word alpha1Scheme("div(phi,alpha1)");
    word alpha1rScheme("div(phirb,alpha1)");

    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(cAlpha*phic, max(phic));
    surfaceScalarField phir(phic*nHatf);

    for (int aCorr=0; aCorr<nAlpha1Corr; aCorr++)
    {
        surfaceScalarField phiAlpha1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alpha1Scheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpha1rScheme),
                alpha1,
                alpha1rScheme
            )
        );

        surfaceScalarField phiAlpha1Rho1
        (
            fvc::flux
            (
                phiAlpha1,
                rho1,
                alpha1Scheme
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
                dimensionSet(1,-3,-1,0,0,0,0),
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
                dimensionSet(1,-3,-1,0,0,0,0),
                0.0
            )
        );

        surfaceScalarField phiRho1 = phi*fvc::interpolate(rho1);

        #ifdef ALPHA1_EXPLICIT_SOLVE
        MULES::explicitSolve
        (
            rho1,
            alpha1,
            phiRho1,
            phiAlpha1Rho1,
            Sp,
            Su,
            1,
            0
        );
        #else
        MULES::implicitSolve
        (
            rho1,
            alpha1,
            phiRho1,
            phiAlpha1Rho1,
            Sp,
            Su,
            1,
            0
        );
        #endif

        volScalarField rRHO1(1.0/rho1);
        surfaceScalarField rRHO1f(fvc::interpolate(rRHO1));

        phiAlpha1 = phiAlpha1Rho1*rRHO1f;
        phiAlpha1Global = phiAlpha1;

        surfaceScalarField rho1f
        (
            "rho1f",
            fvc::interpolate(rho1)
        );

        rhoPhi = phiAlpha1*(rho1f - rho2) + phi*rho2;
    }

    Info<< "alpha1 volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
