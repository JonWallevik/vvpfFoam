// ----------------------------------------------------------------------------
// Conditional compilation with macro definitions (preprocessor directives):
// ----------------------------------------------------------------------------
// Used in gravitySegregation.H:
// To activate the gravity induced segregation of the dispersed phase (alphaD),
// then uncomment the line below:
#define GRAVITY_SEGREGATION
// ----------------------------------------------------------------------------
// Used in particleMigration.H:
// To activate the shear (rate) induced particle migration of the dispersed
// phase (alphaD), then uncomment the line below:
#define PARTICLE_MIGRATION
// ----------------------------------------------------------------------------
// Used in vvpfFoam.H:
// If error analysis is needed, then comment the line below:
// #define ERROR_ANALYSIS
// ----------------------------------------------------------------------------
// Used in vvpfFoam.H, UEqn.H and transportProperties.H:
// To solve for a single reference frame system (e.g. for rheometer with a
// rotating inner geometry), then uncomment the line below:
// #define SINGLE_REFERENCE_FRAME
// ----------------------------------------------------------------------------
// Used in alpha1EqnSubCycle.H:
// To apply alpha1EqnRho.H (in alphaDEqnSubCycle.H) then uncomment, otherwise
// to use alpha1Eqn.H, then comment the line below:
#define ALPHA1RHO_SOLVE
// The default setup is a defined ALPHA1RHO_SOLVE !!!
// ----------------------------------------------------------------------------
// Used in alpha1EqnRho.H:
// To solve the alpha1-equation with MULES::implicitSolve, then comment the
// line below:
#define ALPHA1_EXPLICIT_SOLVE
// The default setup is a defined ALPHA1_EXPLICIT_SOLVE !!!
// ----------------------------------------------------------------------------
// Used in alphaDEqn.H:
// Choose only ONE of the three definitions below:
// - - - - - - - - - - - - - - -
// #define ALPHAD_EXPLICIT_SOLVE
#define ALPHAD_IMPLICIT_SOLVE
// #define ALPHAD_MATRIX_SOLVE
// - - - - - - - - - - - - - - -
// ALPHAD_EXPLICIT_SOLVE = alphaD solved by MULES implicit (FCT)
// ALPHAD_IMPLICIT_SOLVE = alphaD solved by MULES explicit (FCT)
// ALPHAD_MATRIX_SOLVE   = alphaD solved by fvScalarMatrix.solve(): It is
// highly non-conservative and use only for experimentation!
// The default setup is a defined ALPHAD_IMPLICIT_SOLVE !!!
// ----------------------------------------------------------------------------
// Used in pEqn.H:
// To use the pressure residue in the pressure equation (i.e. basically solve
// div(U) = div(R), instead of div(U) = 0), then uncomment the line below:
// #define USE_PRESSURE_RESIDUE
// The default setup is a not-defined USE_PRESSURE_RESIDUE !!!
// ----------------------------------------------------------------------------
