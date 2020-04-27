/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Author
    Dr.ing. Jon Elvar Wallevik, Innovation Center Iceland, Reykjavik, Iceland.

    This work was supported by the Icelandic Research Fund (IRF) --
    grant number 163382-05

Application
    vvpfFoam

Description
    Solver that combines the volume of fluid method (VOF) with with the drift
    flux model (DFM). The mixture is assumed high viscous.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interpolationTable.H"
#include "pimpleControl.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "macroDefinitions.H"
#include "createFunctions.H"
#include "apparentViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "transportProperties.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    scalar sumMagWeightedAverageDivUdt = 0;
    scalar sumWeightedAverageDivUdt = 0;

    scalar sumMagDivRwADt = 0;
    scalar sumDivRwADt = 0;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "rho1MaxMinFields.H"
    rho1 == min(max(rho1,rho1MIN), rho1MAX);

    Info<< "                                                      " << endl;
    Info<< " -----------------------------------------------------" << endl;
    Info<< " Physical start parameters used:                      " << endl;
    Info<< " -----------------------------------------------------" << endl;
    Info<< " Phase 1:                                             " << endl;
    Info<< " mu             " << mu                    << "       " << endl;
    Info<< " tau0           " << tau0                  << "       " << endl;
    Info<< " delta          " << delta                 << "       " << endl;
    Info<< " rho1MAX        " << max(rho1MAX).value()  << " kg/m3 " << endl;
    Info<< " rho1MIN        " << min(rho1MIN).value()  << " kg/m3 " << endl;
    Info<< " - - - - - - - - - - - - - - - - - - - - - - - - - -  " << endl;
    Info<< " rho1 max value " << max(rho1).value()     << " kg/m3 " << endl;
    Info<< " rho1 min value " << min(rho1).value()     << " kg/m3 " << endl;
    Info<< " -----------------------------------------------------" << endl;
    Info<< " Phase 2:                                             " << endl;
    Info<< " eta2           " << eta2                  << "       " << endl;
    Info<< " rho2           " << rho2                  << "       " << endl;
    Info<< " ---------------------------------------------------- " << endl;

    Info<< " Compile options (conditional compilation):           " << endl;

    #ifdef ALPHA1RHO_SOLVE
    Info<< " (1) Using alpha1EqnRho.H";
    #else
    Info<< " (1) Using alpha1Eqn.H";
    #endif

    #ifdef ALPHA1_EXPLICIT_SOLVE
    Info<< " with MULES::explicitSolve " << endl;
    #else
    Info<< " with MULES::implicitSolve " << endl;
    #endif

    #if defined(ALPHAD_EXPLICIT_SOLVE)
    Info<< " (2) Using alphaDEqn.H with MULES::explicitSolve      " << endl;
    #elif defined(ALPHAD_IMPLICIT_SOLVE)
    Info<< " (2) Using alphaDEqn.H with MULES::implicitSolve      " << endl;
    #elif defined(ALPHAD_MATRIX_SOLVE)
    Info<< " (2) Using alphaDEqn.H with fvScalarMatrix            " << endl;
    #else
    Info<< " No solver is set for alphaD in macroDefinitions.H    " << endl;
    abort();
    #endif

    #ifdef GRAVITY_SEGREGATION
    Info<< " (3) Gravity induced segregation (settling) is ON     " << endl;
    #else
    Info<< " (3) Gravity induced segregation (settling) is OFF    " << endl;
    #endif

    #ifdef PARTICLE_MIGRATION
    Info<< " (4) Shear (rate) induced particle migration is ON    " << endl;
    #else
    Info<< " (4) Shear (rate) induced particle migration is OFF   " << endl;
    #endif

    #ifdef SINGLE_REFERENCE_FRAME
    Info<< " (5) Single reference frame (SRF) is used             " << endl;
    Info<< "     Omega = " << omega.value() << " rad/s            " << endl;
    #else
    Info<< " (5) Inertial reference frame is used                 " << endl;
    #endif

    #ifdef ERROR_ANALYSIS
    Info<< " (6) Error analysis is enabled                        " << endl;
    #else
    Info<< " (6) Error analysis is disabled                       " << endl;
    #endif

    Info<< " ---------------------------------------------------- " << endl;
    Info<< " Patches for the current case:                        " << endl;

    {
        const fvPatchList& patchList = mesh.boundary();

        forAll(patchList, patchi)
        {
            const fvPatch& currentPatch = patchList[patchi];
            Info<< "     " << currentPatch.name() << endl;
        }
    }

    Info<< " ---------------------------------------------------- " << endl;
    Info<< "\nStarting time loop\n" << endl;


    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alpha1CourantNo.H"
        #include "alphaDCourantNo.H"
        #include "setDeltaT.H"


        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;


        shearRate = sqrtOfTwo*mag(symm(fvc::grad(U)));
        #include "correctViscosity.H"

        rho1 == (mag(alpha1 - alphaD)*rhoC + alphaD*rhoD)/(alpha1 + delta1);
        rho1 == min(max(rho1,rho1MIN), rho1MAX);

        #include "alpha1EqnSubCycle.H"
        #include "alpha1Interface.H"
        volScalarField sigmaK(sigma*K);

        #include "driftVelocity.H"
        #include "alphaDEqnSubCycle.H"
        #include "enableFieldControl.H"

        while (pimple.loop())
        {
            #include "UEqn.H"

            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            shearRate = sqrtOfTwo*mag(symm(fvc::grad(U)));
            #include "correctViscosity.H"
        }


        #include "enableFieldControl.H"


        #ifdef ERROR_ANALYSIS
        {
            #include "comprsblContErrs.H"
            #include "pEqnResidueErrs.H"
            #include "densityContErrs.H"
        }
        #endif


        Info<< "Averaged shearRateAlpha1 = "
            << shearRateAlpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(shearRateAlpha1) = " << min(shearRateAlpha1).value()
            << "  Max(shearRateAlpha1) = " << max(shearRateAlpha1).value()
            << endl;

        Info<< "Averaged betaD = "
            << betaD.weightedAverage(mesh.Vsc()).value()
            << "  Min(betaD) = " << min(betaD).value()
            << "  Max(betaD) = " << max(betaD).value()
            << endl;

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
