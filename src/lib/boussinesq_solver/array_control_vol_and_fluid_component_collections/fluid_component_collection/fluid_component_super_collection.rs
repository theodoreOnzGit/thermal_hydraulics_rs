// This library was developed for use in my PhD thesis under supervision 
// of Professor Per F. Peterson. It is part of a thermal hydraulics
// library in Rust that is released under the GNU General Public License
// v 3.0. This is partly due to the fact that some of the libraries 
// inherit from GeN-Foam and OpenFOAM, both licensed under GNU General
// Public License v3.0.
//
// As such, the entire library is released under GNU GPL v3.0. It is a strong 
// copyleft license which means you cannot use it in proprietary software.
//
//
/// License
///    This file is part of thermal_hydraulics_rs, a partial library of the
///    thermal hydraulics library written in rust meant to help with the
///    fluid mechanics and heat transfer aspects of the calculations
///     
///    Copyright (C) 2022-2023  Theodore Kay Chen Ong, Singapore Nuclear
///    Research and Safety Initiative, Per F. Peterson, University of 
///    California, Berkeley Thermal Hydraulics Laboratory
///
///    thermal_hydraulics_rs is free software; you can redistribute it and/or modify it
///    under the terms of the GNU General Public License as published by the
///    Free Software Foundation; either version 2 of the License, or (at your
///    option) any later version.
///
///    thermal_hydraulics_rs is distributed in the hope that it will be useful, but WITHOUT
///    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
///    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
///    for more details.
///
///    This library is part of a thermal hydraulics library in rust
///    and contains some code copied from GeN-Foam, and OpenFOAM derivative.
///    This offering is not approved or endorsed by the OpenFOAM Foundation nor
///    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
///    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.
///    Nor is it endorsed by the authors and owners of GeN-Foam.
///
///    You should have received a copy of the GNU General Public License
///    along with this program.  If not, see <http://www.gnu.org/licenses/>.
///
/// © All rights reserved. Theodore Kay Chen Ong,
/// Singapore Nuclear Research and Safety Initiative,
/// Per F. Peterson,
/// University of California, Berkeley Thermal Hydraulics Laboratory
///
/// Main author of the code: Theodore Kay Chen Ong, supervised by
/// Professor Per F. Peterson

use super::fluid_component_collection::*;

/// A struct containing a vector of fluid component collections
#[derive(Debug,Clone,PartialEq)]
pub struct FluidComponentSuperCollection {
    /// this vector contains a collection of fluid component collections
    /// usually, these are in series
    pub fluid_component_super_vector: Vec<FluidComponentCollection>,
}


impl FluidComponentSuperCollection {

    /// returns a copy of the fluid component collection vector
    /// containing immutable elements
    ///
    /// you'll probably need some legwork to create a fresh
    /// object
    ///
    /// the trait object which is given is the FluidComponentCollectionMethods
    /// trait objects
    ///
    /// as such the fluid component super collection may have 
    /// fluid component collections 
    /// and even super collections
    /// in series or parallel or whatever arrangement is desired
    /// as long as it fulfils the fluid component collection methods
    /// trait
    ///
    /// even a single fluid component can behave like a fluid component
    /// collection of 1 item if it fulfils this trait
    ///
    pub fn get_immutable_vector(&self) 
        -> Vec<FluidComponentCollection>{
            self.fluid_component_super_vector.clone()
    }

    /// sets the fluid component collection vector to a specific value
    pub fn set_vector(
        &mut self,
        fluid_component_super_vector: 
        Vec<FluidComponentCollection>){
        self.fluid_component_super_vector = fluid_component_super_vector;
    }


    /// adds a fluid component collection to the super collection

    pub fn add_collection_to_vector(
        &mut self,
        fluid_component_super_vector: Vec<FluidComponentCollection>,
        fluid_component_vector: FluidComponentCollection){

        // first i make a mutable version of the fluid component super vector
        let mut fluid_component_super_vector_mutable =
            fluid_component_super_vector;

        // then i push the pointer to this mutable copy
        fluid_component_super_vector_mutable.push(fluid_component_vector);

        // next i set the fluid component vector
        self.set_vector(fluid_component_super_vector_mutable);

    }

    /// removes a fluid component collection by index from the super collection

    pub fn remove_collection_by_index(&mut self,
              fluid_component_super_vector: 
              Vec<FluidComponentCollection>,
              component_index: usize){

        // first i make a mutable copy of the component vector
        let mut fluid_component_super_vector_mutable =
            fluid_component_super_vector;

        // i remove the index from the vector 
        // (note that there may be a case where the vector is smaller than
        // the given index),
        // however, the remove method already has a panic if the 
        // vector is shorter than the given index

        fluid_component_super_vector_mutable.remove(component_index);

        // next i set the fluid component vector
        self.set_vector(fluid_component_super_vector_mutable);
    }

    /// returns read only a pointer of the fluid component collection
    /// given an index

    pub fn get_collection_by_index(
        &mut self,
        component_index: usize) -> FluidComponentCollection{

        // first let's access the fluid component super vector

        let fluid_component_super_vector =
            self.get_immutable_vector();

        // then i access a particular super collection

        let fluid_component_collection_pointer = 
            fluid_component_super_vector[component_index].clone();

        return fluid_component_collection_pointer;

    }


    /// updates the fluid component collection at the specified
    /// index with a fluid component collection supplied by the user

    pub fn update_collection_by_index(
        &mut self,
        component_index: usize,
        fluid_component_super_vector: Vec<FluidComponentCollection>,
        fluid_component_collection: FluidComponentCollection){

        // first i make a mutable copy of the component vector
        let mut fluid_component_super_vector_mutable =
            fluid_component_super_vector;

        // then i change the pointer in this mutable copy
        fluid_component_super_vector_mutable[component_index]
            = fluid_component_collection;

        // next i set the fluid component vector
        self.set_vector(fluid_component_super_vector_mutable);
    }



}

/// the default is to provide an empty vector
impl Default for FluidComponentSuperCollection {
    fn default() -> Self {
        Self { fluid_component_super_vector: vec![] }
    }
}
