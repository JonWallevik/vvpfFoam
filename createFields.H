    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha1\n" << endl;
    volScalarField alpha1
    (
        IOobject
        (
            "alpha1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField shearRate
    (
        IOobject
        (
            "shearRate",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "shearRate",
            dimensionSet(0,0,-1,0,0,0,0),
            scalar(0)
        )
    );

    volScalarField shearRateAlpha1
    (
        IOobject
        (
            "shearRateAlpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        shearRate
    );

    volScalarField etaEff
    (
        IOobject
        (
            "etaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "etaEff",
            dimensionSet(1,-1,-1,0,0,0,0),
            scalar(0)
        )
    );

    volVectorField gradEtaEff
    (
        IOobject
        (
            "gradEtaEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "gradEtaEff",
            dimensionSet(1,-2,-1,0,0,0,0),
            vector::zero
        )
    );

    Info<< "Reading field alphaD\n" << endl;
    volScalarField alphaD
    (
        IOobject
        (
            "alphaD",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField betaD
    (
        IOobject
        (
            "betaD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        alphaD/(alpha1 + delta1)
    );

    volScalarField wallBetaD
    (
        IOobject
        (
            "wallBetaD",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        betaD
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    volVectorField Ualpha1
    (
        IOobject
        (
            "Ualpha1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "Ualpha1",
            dimensionSet(0,1,-1,0,0,0,0),
            vector::zero
        )
    );

    // For single reference frame calculations (inertial velocity):
    #ifdef SINGLE_REFERENCE_FRAME
    volVectorField Uin
    (
        IOobject
        (
            "Uin",
             runTime.timeName(),
             mesh,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "Uin",
            U.dimensions(),
            vector::zero
        )
    );
    #endif

    volVectorField VdjGR
    (
        IOobject
        (
            "VdjGR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "VdjGR",
            dimensionSet(0,1,-1,0,0,0,0),
            vector::zero
        )
    );

    surfaceScalarField phiVdjGR
    (
        IOobject
        (
            "phiVdjGR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(VdjGR) & mesh.Sf()
    );

    volVectorField VdjSR
    (
        IOobject
        (
            "VdjSR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "VdjSR",
            dimensionSet(0,1,-1,0,0,0,0),
            vector::zero
        )
    );

    surfaceScalarField phiVdjSR
    (
        IOobject
        (
            "phiVdjSR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(VdjSR) & mesh.Sf()
    );

    volVectorField Vdj
    (
        IOobject
        (
            "Vdj",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        VdjGR + VdjSR
    );

    surfaceScalarField phiVdj
    (
        IOobject
        (
            "phiVdj",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(VdjGR + VdjSR) & mesh.Sf()
    );

    volScalarField rho1
    (
        IOobject
        (
            "rho1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (scalar(1) - alphaD)*rhoC + alphaD*rhoD
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1*rho1 + (scalar(1) - alpha1)*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();

    surfaceScalarField phiAlpha1Global
    (
        IOobject
        (
            "alpha1*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(alpha1)*phi
    );

    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );

    surfaceScalarField rhoPhiVdjGR
    (
        IOobject
        (
            "rhoPhiVdjGR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phiVdjGR
    );

    surfaceScalarField rhoPhiVdjSR
    (
        IOobject
        (
            "rhoPhiVdjSR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phiVdjSR
    );

    surfaceScalarField rhoPhiVdj
    (
        IOobject
        (
            "rhoPhiVdj",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*(phiVdjGR + phiVdjSR)
    );

    surfaceScalarField rhoPhiAlphaD
    (
        IOobject
        (
            "rho1*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoPhi
    );

    #include "readGravitationalAcceleration.H"

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    scalar cAlpha
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cAlpha")
        )
    );

    surfaceScalarField nHatf
    (
        IOobject
        (
            "nHatf",
            alpha1.time().timeName(),
            alpha1.mesh()
        ),
        alpha1.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    );

    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    );

    volScalarField K
    (
        IOobject
        (
            "K",
            alpha1.time().timeName(),
            alpha1.mesh()
        ),
        alpha1.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    );

    // fv::IOoptionList fvOptions(mesh);
