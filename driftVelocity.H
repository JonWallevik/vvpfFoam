{
    // lower and upper bounds:
    const dimensionedScalar lowerCrit("lowerCrit", dimless, scalar(0.02));
    const dimensionedScalar upperCrit("upperCrit", dimless, scalar(0.25));

    // Note, criteriaD has been determined when alphaDMIN = 0. A different
    // setup might be needed for criteriaD, if alphaDMIN > 0:
    const dimensionedScalar criteriaD
    (
        0.100*(alphaDMAX - alphaDMIN) + alphaDMIN
    );

    // zero velocity:
    const dimensionedVector zeroVelocity
    (
        "zeroVelocity",
        dimensionSet(0,1,-1,0,0,0,0),
        vector::zero
    );

    // see alpha1Interface.H:
    volVectorField gradAlpha1(fvc::grad(alpha1));
    surfaceVectorField gradAlpha1f(fvc::interpolate(gradAlpha1));
    surfaceVectorField interfaceNormal(gradAlpha1f/(mag(gradAlpha1f) + deltaN));

    #include "gravitySegregation.H"
    #include "particleMigration.H"

    phiVdj = phiVdjGR + phiVdjSR;
    Vdj    = VdjGR    + VdjSR;
}
