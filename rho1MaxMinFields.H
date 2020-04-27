const dimensionedScalar rho1MAXscalar =
    (scalar(1) - alphaDMAX)*rhoC + alphaDMAX*rhoD;

const dimensionedScalar rho1MINscalar =
    (scalar(1) - alphaDMIN)*rhoC + alphaDMIN*rhoD;


const volScalarField rho1MAX
(
    IOobject
    (
        "rho1MAX",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "rho1MAX",
        dimensionSet(1,-3,0,0,0,0,0),
        rho1MAXscalar.value()
    )
);


const volScalarField rho1MIN
(
    IOobject
    (
        "rho1MIN",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "rho1MIN",
        dimensionSet(1,-3,0,0,0,0,0),
        rho1MINscalar.value()
    )
);
