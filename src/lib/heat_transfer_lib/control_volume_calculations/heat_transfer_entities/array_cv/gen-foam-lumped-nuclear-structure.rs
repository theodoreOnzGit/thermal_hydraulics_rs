// I'm translating lumpedNuclearStructure.C from GeN-Foam into Rust
/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2022 OpenCFD Ltd.         |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/
//
//#include "lumpedNuclearStructure.H"
//#include "structure.H"
//#include "addToRunTimeSelectionTable.H"
//#include "SquareMatrix.H"
//
//// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
//
//namespace Foam
//{
//namespace powerModels
//{
//    defineTypeNameAndDebug(lumpedNuclearStructure, 0);
//    addToRunTimeSelectionTable
//    (
//        powerModel, 
//        lumpedNuclearStructure, 
//        powerModels
//    );
//}
//}
//
//
//// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//
//Foam::powerModels::lumpedNuclearStructure::lumpedNuclearStructure
//(
//    structure& structureRef,
//    const dictionary& dicts
//)
//:
//    powerModel
//    (
//        structureRef,
//        dicts
//    ),
//    T_
//    (
//        IOobject
//        (
//            "T."+typeName,
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh_.cells().size()
//    ),
//    Tmax_
//    (
//        IOobject
//        (
//            "Tmax."+typeName,
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar("", dimTemperature, 0),
//        zeroGradientFvPatchScalarField::typeName
//    ),
//    Tsurface_
//    (
//        IOobject
//        (
//            "Tsurface."+typeName,
//            mesh_.time().timeName(),
//            mesh_,
//            IOobject::READ_IF_PRESENT,
//            IOobject::AUTO_WRITE
//        ),
//        mesh_,
//        dimensionedScalar("", dimTemperature, 0),
//        zeroGradientFvPatchScalarField::typeName
//    ),
//    fractionOfPowerFromNeutronics_(0),
//    nodesNumber_(0),
//    nodeFuel_(0),
//    nodeClad_(0),
//    Hs_(0),
//    rhoCp_(0),
//    volFraction_(0),
//    qFraction_(0),
//    cellToRegion_(mesh_.cells().size(), 0),
//    regionIndexToRegionName_(0)
//{   
//    this->setInterfacialArea();
//    structure_.setRegionField(*this, structureRef.powerDensityNeutronics(), "powerDensity");
//
//    typedef IOFieldField<Field, scalar> scalarFieldField;
//    bool foundT
//    (
//        T_.typeHeaderOk<scalarFieldField>(true)
//    );
//
//    scalarList T0(0);
//    forAll(this->toc(), regioni)
//    {
//        word region(this->toc()[regioni]);
//        const dictionary& dict(this->subDict(region));       
//        //- Setup cellToRegion_ mapping
//        const labelList& regionCells
//        (
//            structure_.cellLists()[region]
//        );
//        forAll(regionCells, i)
//        {
//            label celli(regionCells[i]);
//            cellToRegion_[celli] = regioni;
//        }
//
//        //- Add to regionIndexToRegionName_ mapping
//        regionIndexToRegionName_.append(region);
//
//        //- Add to regionIndexToRegionName_ mapping
//        regionIndexToRegionName_.append(region);
//
//        //- Read region dict entries
//        scalar fractionOfPowerFromNeutronics(dict.lookupOrDefault<scalar>("fractionOfPowerFromNeutronics",1.0));
//        label nodesNumber(dict.get<label>("nodesNumber"));
//        label nodeFuel(dict.get<label>("nodeFuel"));
//        label nodeClad(dict.get<label>("nodeClad"));
//        scalarList Hs(dict.get<scalarList>("heatConductances"));
//        scalarList rhoCp(dict.get<scalarList>("rhoCp"));
//        scalarList volFraction(dict.get<scalarList>("volumeFractions"));
//        scalarList qFraction(dict.get<scalarList>("powerFractions"));
//
//        if (!foundT)
//        {
//            T0.append(dict.get<scalar>("T0"));
//        }
//
//        //- Fill in lists for this region
//        fractionOfPowerFromNeutronics_.append(fractionOfPowerFromNeutronics), 
//        nodesNumber_.append(nodesNumber),
//        nodeFuel_.append(nodeFuel),
//        nodeClad_.append(nodeClad),
//        Hs_.append(Hs);
//        rhoCp_.append(rhoCp);
//        volFraction_.append(volFraction);
//        qFraction_.append(qFraction);
//    }
//    //- If T not found, init it from dictionary value T0 read previously
//    if (!foundT)
//    {
//        forAll(mesh_.cells(), i)
//        {
//            T_.set(i, new Field<scalar>(0, 0));
//        }
//
//        Info<< "Reading lumpedNuclearStructure initial temperatures from "
//            << "dictionary" << endl;
//
//        forAll(this->cellList_, i)
//        {
//            label celli(this->cellList_[i]);
//            label regioni(cellToRegion_[celli]);
//            T_.set(celli, new Field<scalar>(nodesNumber_[regioni], 0));
//            forAll(T_[celli], subCelli)
//            {
//                T_[celli][subCelli] = T0[regioni];
//            }
//        }
//    }
//    else
//    {
//        Info<< "Setting lumpedNuclearStructure initial temperatures from "
//                << T_.name() << endl;
//    }
//    
//    //- Set I/O fields and compute initial scalar max, min
//    forAll(this->cellList_, i)
//    {
//        label celli(this->cellList_[i]);
//        label regioni(cellToRegion_[celli]);
//        const scalarField& T(T_[celli]);
//
//        // Update average fuel and clad temp used for coupling
//        this->structureRef().TFuelAv()[celli] = T[nodeFuel_[regioni]];
//        this->structureRef().TCladAv()[celli] = T[nodeClad_[regioni]];
//    }
//}
//
//
//// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
//Foam::powerModels::lumpedNuclearStructure::~lumpedNuclearStructure()
//{}
//
//// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//
//void 
//Foam::powerModels::lumpedNuclearStructure::updateLocalTemperatureProfile
//(
//    const label& celli,
//    const scalar& HTSumi,
//    const scalar& HSumi
//)
//{
//    //-
//    scalarField& T(T_[celli]);
//
//    //- Read region values
//    const label& regioni(cellToRegion_[celli]);
//    const word& region(regionIndexToRegionName_[regioni]);
//    const scalar& fractionOfPowerFromNeutronics(fractionOfPowerFromNeutronics_[regioni]);
//    const label& nodesNumber(nodesNumber_[regioni]);
//    const label& nodeFuel(nodeFuel_[regioni]);
//    const label& nodeClad(nodeClad_[regioni]);
//    const scalarList& Hs(Hs_[regioni]);
//    const scalarList& rhoCp(rhoCp_[regioni]);
//    const scalarList& volFraction(volFraction_[regioni]);
//    const scalarList& qFraction(qFraction_[regioni]);
//    
//    const scalarField& TOld = T_.oldTime()[celli];
//
//    const scalar& iA = this->iA_[celli];
//
//    //- Update power density
//    const scalar& qRef(structure_.powerDensityNeutronics()[celli]);
//    scalar q = qRef * fractionOfPowerFromNeutronics;
//
//    //- Recurrent quantities
//    scalar dt(mesh_.time().deltaT().value());
//    scalar Tcool(HTSumi / max(HSumi,SMALL));
//    scalar Hcool(HSumi*iA/this->alpha_[celli]);
//
//    if(nodesNumber>1)
//    {
//        //- Init matrix, source
//        SquareMatrix<scalar> M(nodesNumber, nodesNumber, Foam::zero());
//        List<scalar> S(nodesNumber, 0.0);
//
//        //- Construct matrix, source
//        {
//            //- Set "zeroGradient" BC at innermost node
//            {
//                M[0][1] =   -Hs[1];
//                M[0][0] =   volFraction[0] * rhoCp[0] / dt + Hs[1];
//                S[0] =      q * qFraction[0]  + TOld[0] * volFraction[0] * rhoCp[0] / dt;
//            }
//
//            // Bulk
//            if(nodesNumber>2)
//            {
//                for (int i = 1; i < nodesNumber-1; i++)
//                {
//                    M[i][i+1] =     -Hs[i+1];
//                    M[i][i-1] =     -Hs[i];
//                    M[i][i] =       volFraction[i] * rhoCp[i] / dt + Hs[i+1] + Hs[i];
//                    S[i] =          q * qFraction[i] + TOld[i] * volFraction[i] * rhoCp[i] / dt;
//                }
//            }
//
//            //- Outer surface, convective BC with fluid(s) wetting the pin
//            {
//                label i(nodesNumber-1);
//                scalar HtoCool(Hs[i+1]*Hcool/(Hs[i+1]+Hcool)); //total H from last node to coolant
//                M[i][i-1] =     -Hs[i];
//                M[i][i] =       volFraction[i] * rhoCp[i] / dt + Hs[i] + HtoCool;
//                S[i] =          q * qFraction[i] 
//                                + TOld[i] * volFraction[i] * rhoCp[i] / dt 
//                                + HtoCool * Tcool;       
//            }
//        }
//
//        //- Solve linear system
//        solve(T, M, S);
//    }
//    else
//    {
//        scalar HtoCool(Hs[1]*Hcool/(Hs[1]+Hcool)); //total H from last node to coolant
//        scalar M(volFraction[0] * rhoCp[0] / dt + HtoCool);
//        scalar S
//                (
//                    q * qFraction[0] 
//                    + TOld[0] * volFraction[0] * rhoCp[0] / dt 
//                    + HtoCool * Tcool
//                );
//        T = S/M;
//    }
//
//    //- Set fields (max and outer)
//    Tmax_[celli] = T[0] + q * qFraction[0] / Hs[0];
//    Tsurface_[celli] = (Hs[nodesNumber] * T[nodesNumber-1] + Hcool * Tcool ) 
//                        / (Hs[nodesNumber] + Hcool);
//
//    /*
//    Info << "T " << T << endl;    
//    Info << "Tmax_[celli] " << Tmax_[celli] << endl;
//    Info << "Tsurface_[celli] " << Tsurface_[celli] << endl;
//    Info << "Tcool " << Tcool << endl;   
//    */
//
//    // Update average fuel and clad temp used for coupling
//    this->structureRef().TFuelAv()[celli] = T[nodeFuel_[regioni]];
//    this->structureRef().TCladAv()[celli] = T[nodeClad_[regioni]];
//}
//
//
//void Foam::powerModels::lumpedNuclearStructure::correct
//(
//    const volScalarField& HTSum,  // == SUM_j [htc_j*T_j*frac_j]
//    const volScalarField& HSum    // == SUM_j [htc_j*frac_j]
//)
//{
//    
//    //- Update temperatures cell-by-cell
//    forAll(this->cellList_, i)
//    {
//        label celli(this->cellList_[i]);
//        updateLocalTemperatureProfile(celli, HTSum[celli], HSum[celli]);
//    }
//
//}
//
//void Foam::powerModels::lumpedNuclearStructure::correctT(volScalarField& T) const
//{
//    //- Set T to  surface temperature
//    forAll(cellList_, i)
//    {
//        label celli(cellList_[i]);
//        T[celli] = Tsurface_[celli];
//    }
//}
//
//// ************************************************************************* //
//

use ndarray::*;
use ndarray_linalg::*;
use uom::si::f64::*;
use uom::si::thermal_conductance::watt_per_kelvin;
use uom::si::thermodynamic_temperature::kelvin;
use uom::si::time::second;
use uom::ConstZero;

// 
#[inline]
#[allow(non_snake_case)]
pub (in crate::heat_transfer_lib::control_volume_calculations)
fn matrix_construction(){
    // translating:
    //    scalar dt(mesh_.time().deltaT().value());
    //    scalar Tcool(HTSumi / max(HSumi,SMALL));
    //    scalar Hcool(HSumi*iA/this->alpha_[celli]);
    //
    //
    // there is no coolant temperature so i'll just call it boundary 
    // condition B
    //
    // also boundary condition B has a thermal conductance
    // which I think should be user determined

    let dt: Time = Time::new::<second>(0.1);
    let boundary_condition_b_temperature: ThermodynamicTemperature = 
    ThermodynamicTemperature::new::<kelvin>(283.0);
    let boundary_condition_b_conductance: ThermalConductance = 
    ThermalConductance::new::<watt_per_kelvin>(0.1);


    // translating: 
    //    if(nodesNumber>1)
    //    
    let nodesNumber: usize = 10;
    if nodesNumber > 1 {
        //    translating:
        //    {
        //        //- Init matrix, source
        //        SquareMatrix<scalar> M(nodesNumber, nodesNumber, Foam::zero());
        //        List<scalar> S(nodesNumber, 0.0);

        let M: Array2<ThermalConductance> = 
        Array::zeros((nodesNumber, nodesNumber));
        let S: Array1<Power> = Array::zeros(nodesNumber);

        // translating:
        //        //- Construct matrix, source
        //        {
        //            //- Set "zeroGradient" BC at innermost node
        //            {
        //                M[0][1] =   -Hs[1];
        //                M[0][0] =   volFraction[0] * rhoCp[0] / dt + Hs[1];
        //                S[0] =      q * qFraction[0]  + TOld[0] * volFraction[0] 
        //                * rhoCp[0] / dt;
        //            }
        //
        // now, Hs is the conductance vector, but in GeN-Foam this is done 
        // on a per unit volume basis
        

    } else {

    }
    //
    //            // Bulk
    //            if(nodesNumber>2)
    //            {
    //                for (int i = 1; i < nodesNumber-1; i++)
    //                {
    //                    M[i][i+1] =     -Hs[i+1];
    //                    M[i][i-1] =     -Hs[i];
    //                    M[i][i] =       volFraction[i] * rhoCp[i] / dt + Hs[i+1] + Hs[i];
    //                    S[i] =          q * qFraction[i] + TOld[i] * volFraction[i] * rhoCp[i] / dt;
    //                }
    //            }
    //
    //            //- Outer surface, convective BC with fluid(s) wetting the pin
    //            {
    //                label i(nodesNumber-1);
    //                scalar HtoCool(Hs[i+1]*Hcool/(Hs[i+1]+Hcool)); //total H from last node to coolant
    //                M[i][i-1] =     -Hs[i];
    //                M[i][i] =       volFraction[i] * rhoCp[i] / dt + Hs[i] + HtoCool;
    //                S[i] =          q * qFraction[i] 
    //                                + TOld[i] * volFraction[i] * rhoCp[i] / dt 
    //                                + HtoCool * Tcool;       
    //            }
    //        }
    //
    //        //- Solve linear system
    //        solve(T, M, S);
    //    }
    //    else
    //    {
    //        scalar HtoCool(Hs[1]*Hcool/(Hs[1]+Hcool)); //total H from last node to coolant
    //        scalar M(volFraction[0] * rhoCp[0] / dt + HtoCool);
    //        scalar S
    //                (
    //                    q * qFraction[0] 
    //                    + TOld[0] * volFraction[0] * rhoCp[0] / dt 
    //                    + HtoCool * Tcool
    //                );
    //        T = S/M;
    //    }

    
}

#[inline]
fn solve_conductance_matrix_power_vector(
    thermal_conductance_matrix: Array2<ThermalConductance>,
    power_vector: Array1<Power>)
-> Result<Array1<ThermodynamicTemperature>, error::LinalgError>{

    // I can of course convert it into f64 types 
    //
    //

    let get_value_conductance = |conductance: &ThermalConductance| {
        return conductance.value;
    };

    let get_value_power = |power: &Power| {
        return power.value;
    };

    // i'm allowing non snake case so that the syntax is the same as 
    // GeN-Foam
    #[allow(non_snake_case)]
    let M: Array2<f64> = 
    thermal_conductance_matrix.map(get_value_conductance);

    #[allow(non_snake_case)]
    let S: Array1<f64> = power_vector.map(get_value_power);

    // now for the raw temperature matrix 

    #[allow(non_snake_case)]
    let T: Array1<f64> = M.solve(&S)?;

    // To check for unit safety, I can just perform one calc

    let _unit_check: Power = 
    power_vector[0] + 
    thermal_conductance_matrix[[0,0]] 
    * ThermodynamicTemperature::ZERO;

    // now map T back to a ThermodynamicTemperature
    // T is already a ThermodynamicTemperature, so don't need to manually 
    // convert, do it in kelvin 
    //

    let convert_f64_to_kelvin_temperature = |float: &f64| {
        return ThermodynamicTemperature::new::<kelvin>(*float);
    };


    // this is the last step
    let temperature_vector: Array1<ThermodynamicTemperature> 
    = T.map(convert_f64_to_kelvin_temperature);

    return Ok(temperature_vector);
}
