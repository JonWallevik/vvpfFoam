const scalar PI = constant::mathematical::pi;
const scalar convertToRad = PI/180.0;

/*---------------------------------------------------------------------------*\
The term "varPhi" (see $\varphi$ in LaTex) designates volume fraction of the
dispersed phase. The term varPhi can either be alphaD or betaD, depending on
user preference. The modification is made in the following source codes:
correctViscosity.H, gravitySegregation.H and particleMigration.H.
Note, making the same choice in all the three codes is not necessary! The best
choice will depend on the specific case.
\*---------------------------------------------------------------------------*/

// A test function:
volScalarField slowDown1
(
    const volScalarField& varPhi,
    const dimensionedScalar& varPhiMIN,
    const dimensionedScalar& varPhiMAX,
    const dimensionedScalar& varPhiNORM
)
{
    return (0.5 + (1.0/PI)*atan(-15.0
         * (mag(varPhi - varPhiMIN)/varPhiMAX
         - 0.5*(varPhiMAX + varPhiNORM)/varPhiMAX)));
}

/*---------------------------------------------------------------------------*\
The function below is used to prevent settlement by gravity into cells that
already have reached the maximum allowed volume fraction, i.e. varPhiMAX.
It is a arctan function, mostly equal to 1, with a fast reduction towards zero
as varPhi approaches varPhiMAX.
octave: varPhiMIN = 0; varPhiMAX = 0.4; varPhi = [varPhiMIN:0.0001:varPhiMAX];
slowDown2 = (2/pi).*atan(-400.0.*((varPhi - varPhiMIN)./varPhiMAX -
            (varPhiMAX - varPhiMIN)/varPhiMAX));
\*---------------------------------------------------------------------------*/

volScalarField slowDown2
(
    const volScalarField& varPhi,
    const dimensionedScalar& varPhiMIN,
    const dimensionedScalar& varPhiMAX
)
{
    return ((2.0/PI)*atan(-400.0*(mag(varPhi - varPhiMIN)/varPhiMAX
         - mag(varPhiMAX - varPhiMIN)/varPhiMAX)));
}

// A test function:
volScalarField slowDown3
(
    const volScalarField& alpha1,
    const volScalarField& varPhi,
    const dimensionedScalar& varPhiMIN,
    const dimensionedScalar& varPhiMAX,
    const dimensionedScalar& varPhiNORM
)
{
    return (1.00*(mag(varPhi - varPhiMIN)/varPhiMAX)
         * mag(pow((1.0 - mag(varPhi - varPhiMIN)/varPhiMAX), 0.5))
         * mag(alpha1));
}

/*---------------------------------------------------------------------------*\
The function below is used to prevent shear (rate) induced particle migration
into cells that have already reached the maximum allowed volume fraction, i.e.
varPhiMAX. That is, the function advances towards zero as varPhi approaches
varPhiMAX. The function also advances towards zero as varPhi approaches
varPhiMIN. More precisely, this is a symmetric function, mostly equal to 1,
with a fast reduction towards zero as varPhi approaches varPhiMIN or varPhiMAX.
octave: varPhiMIN = 0; varPhiMAX = 0.4; varPhi = [varPhiMIN:0.0001:varPhiMAX];
slowDown4 = (2/pi).*atan(-400.0.*((varPhi - varPhiMIN)./varPhiMAX - ...
            (varPhiMAX - varPhiMIN)/varPhiMAX))* ...
            (2/pi).*atan(400.0.*((varPhi - varPhiMIN)/varPhiMAX));
\*---------------------------------------------------------------------------*/

volScalarField slowDown4
(
    const volScalarField& varPhi,
    const dimensionedScalar& varPhiMIN,
    const dimensionedScalar& varPhiMAX
)
{
    return ((2.0/PI)*atan(-400.0*(mag(varPhi - varPhiMIN)/varPhiMAX
         - mag(varPhiMAX - varPhiMIN)/varPhiMAX))
         * (2.0/PI)*atan(400.0*(mag(varPhi - varPhiMIN)/varPhiMAX)));
}
