{
    const fvMesh& mesh = alpha1.mesh();
    const surfaceVectorField& Sf = mesh.Sf();

    // Cell gradient of alpha1
    const volVectorField gradAlpha1(fvc::grad(alpha1));

    // Interpolated face-gradient of alpha1
    surfaceVectorField gradAlpha1f(fvc::interpolate(gradAlpha1));

    // Face unit interface normal
    surfaceVectorField nHatfv(gradAlpha1f/(mag(gradAlpha1f) + deltaN));

    surfaceVectorField::GeometricBoundaryField nHatb =
        nHatfv.boundaryField();
    surfaceVectorField::GeometricBoundaryField gradAlpha1fb =
        gradAlpha1f.boundaryField();

    const volScalarField::GeometricBoundaryField& abf = alpha1.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();


    forAll(boundary, patchi)
    {
        /*
        Info<< "isA<alphaContactAngleFvPatchScalarField>(abf[patchi])"
            << isA<alphaContactAngleFvPatchScalarField>(abf[patchi]) << endl;

        const fvPatch& currentPatch = boundary[patchi];
        Info<< "     " << currentPatch.name() << endl;
        */

        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(Foam::cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = Foam::cos(Foam::acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlpha1fb[patchi]);
            acap.evaluate();
       }
    }

    nHatfv.boundaryField() = nHatb;
    gradAlpha1f.boundaryField() = gradAlpha1fb;

    // Face unit interface normal flux
    nHatf = nHatfv & Sf;

    // Simple expression for curvature
    K = -fvc::div(nHatf);
}
