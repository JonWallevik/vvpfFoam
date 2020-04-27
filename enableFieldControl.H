gradEtaEff = fvc::grad(etaEff);

forAll(alpha1.internalField(), celli)
{
    if (alpha1[celli] < alpha1EtaCriteria.value())
    {
        gradEtaEff[celli] = vector::zero;
    }
}


shearRateAlpha1 == alpha1*shearRate; // smoother!
Ualpha1 == alpha1*U;

betaD == alphaD/(alpha1 + delta1);
betaD == min(max(betaD, alphaDMIN), alphaDMAX);

wallBetaD.internalField() = betaD.internalField();


{
    const fvPatchList& patchList = mesh.boundary();

    forAll(patchList, patchi)
    {
        const fvPatch& currentPatch = patchList[patchi];

        forAll(currentPatch, facei)
        {
            label faceCelli = currentPatch.faceCells()[facei];
            wallBetaD.boundaryField()[patchi][facei] = wallBetaD[faceCelli];
        }
    }
}


// For single reference frame calculations (inertial velocity):
#ifdef SINGLE_REFERENCE_FRAME
Uin = U + (Omega ^ mesh.C());
// Uin = alpha1*(U + (Omega ^ mesh.C()));
#endif


forAll(alpha1.internalField(), celli)
{
    if (alpha1[celli] < alpha1Criteria.value())
    {
        shearRateAlpha1[celli] = 0.0; // cut off!
        betaD[celli] = 0.0;
        wallBetaD[celli] = 0.0;
        rho1[celli] = deltaRho1.value();
        Ualpha1[celli] = vector::zero;
        #ifdef SINGLE_REFERENCE_FRAME
        Uin[celli] = vector::zero;
        #endif
    }
}


forAll(mesh.boundary(), patchi)
{
    forAll(alpha1.boundaryField()[patchi], facei)
    {
        if (alpha1.boundaryField()[patchi][facei] < alpha1Criteria.value())
        {
            shearRateAlpha1.boundaryField()[patchi][facei] = 0.0;
            betaD.boundaryField()[patchi][facei] = 0.0;
            wallBetaD.boundaryField()[patchi][facei] = 0.0;
            rho1.boundaryField()[patchi][facei] = deltaRho1.value();
            Ualpha1.boundaryField()[patchi][facei] = vector::zero;
            #ifdef SINGLE_REFERENCE_FRAME
            Uin.boundaryField()[patchi][facei] = vector::zero;
            #endif
        }
    }
}
