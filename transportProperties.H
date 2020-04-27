// ------------------------------------------------------------------------- //

const dimensionedScalar sqrtOfTwo("sqrtOfTwo", dimless, std::sqrt(2.0));

// Used in enableFieldControl.H as a search criteria. See enableFieldControl.H:
const dimensionedScalar alpha1Criteria("alpha1Criteria", dimless, scalar(0.35));

const dimensionedScalar alpha1EtaCriteria
(
    "alpha1EtaCriteria",
    dimless,
    scalar(0.1)
);

// Used for "betaD == alphaD/(alpha1 + delta1)", c.f. enableFieldControl.H:
const dimensionedScalar delta1("delta1", dimless, scalar(1.0e-4));

// Used in enableFieldControl.H, to avoid singularities in "/rho1"-calculations.
// For example, if mixture density is > 2000 kg/m3, the value of deltaRho1 can
// be as high as 1.0 kg/m3. Below 1.0e-3 is not recommended:
const dimensionedScalar deltaRho1("deltaRho1", dimDensity, scalar(1.0e-2));

// Used in particleMigration.H, for gradEtaEff/(etaEff + deltaEta):
const dimensionedScalar deltaEta
(
    "deltaEta",
    dimensionSet(1,-1,-1,0,0,0,0),
    scalar(1.0e-2)
);

// ------------------------------------------------------------------------- //
IOdictionary physicalProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// ------------- Rotational properties -------------

// For single reference frame calculations:
#ifdef SINGLE_REFERENCE_FRAME
dimensionedScalar omega
(
    physicalProperties.lookup("omega")
);

// Omega is in [rad/s]:
dimensionedVector Omega
(
    "Omega",
    dimensionSet(0, 0, -1, 0, 0),
    vector(0, 0, omega.value())
    //    (X, Y, Z)
);
#endif

// ------------- Phase 1 properties -------------
// The following is an example of interaction with case_dir/constant/
// transportProperties and is relative to "return viscous_5();" in
// in the source code apparentViscosity.H.

// (Plastic) viscosity of phase 1:
dimensionedScalar mu
(
    // physicalProperties.lookup("mu")
    physicalProperties.subDict("PHASE1").subDict("RHEOLOGY").lookup("mu")
);

// Yield stress of phase 1:
dimensionedScalar tau0
(
    // physicalProperties.lookup("tau0")
    physicalProperties.subDict("PHASE1").subDict("RHEOLOGY").lookup("tau0")
);

// The regularization parameter:
dimensionedScalar delta
(
    // physicalProperties.lookup("delta")
    physicalProperties.subDict("PHASE1").subDict("RHEOLOGY").lookup("delta")
);

// Density of the continuous phase:
dimensionedScalar rhoC
(
    physicalProperties.subDict("PHASE1").lookup("rhoContinous")
);

// Density of the dispersed phase:
dimensionedScalar rhoD
(
    physicalProperties.subDict("PHASE1").lookup("rhoDispersed")
);

// Characteristic particle diameter of the dispersed phase:
dimensionedScalar Da
(
    physicalProperties.subDict("PHASE1").lookup("Da")
);

// Max consentration allowed for the dispersed phase:
dimensionedScalar alphaDMAX
(
    physicalProperties.subDict("PHASE1").lookup("alphaDMAX")
);

// Min consentration allowed for the dispersed phase:
dimensionedScalar alphaDMIN
(
    physicalProperties.subDict("PHASE1").lookup("alphaDMIN")
);

// ------------- Phase 2 properties -------------

// Apparent viscosity of phase 2 (see apparentViscosity.H about usage):
dimensionedScalar eta2
(
    physicalProperties.subDict("PHASE2").lookup("mu")
);

// Density of phase 2:
dimensionedScalar rho2
(
    physicalProperties.subDict("PHASE2").lookup("rho")
);

// ----------- Phases 1-2 properties ------------

// Surface tension of the interface:
dimensionedScalar sigma
(
    physicalProperties.lookup("sigma")
);

// ------------------------------------------------
