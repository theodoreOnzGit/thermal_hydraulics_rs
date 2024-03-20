use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_solid_array_with_lateral_coupling::SolidColumn;
use crate::boussinesq_solver::boussinesq_thermophysical_properties::Material;
use crate::boussinesq_solver::single_control_vol::SingleCVNode;
use crate::boussinesq_solver::array_control_vol_and_fluid_component_collections::one_d_fluid_array_with_lateral_coupling::FluidArray;
use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use uom::si::f64::*;


/// Contains Types of Control Volumes (CVs)
#[derive(Debug,Clone,PartialEq)]
pub enum CVType {
    /// This CV is the most basic,  it can be represented by a single 
    /// point or node
    SingleCV(SingleCVNode),
    /// Array CVs are collections of SingleCVs, 
    /// or discretised arrays of control volumes with SingleCVNodes 
    /// attached to either end
    /// but do not require the 
    /// user to manually specify the connections between the SingleCVs
    /// This is for fluid arrays, where there is advection through 
    /// the array
    FluidArrayCV(FluidArray),
    /// Array CVs are collections of SingleCVs, 
    /// or discretised arrays of control volumes with SingleCVNodes 
    /// attached to either end
    /// but do not require the 
    /// user to manually specify the connections between the SingleCVs
    /// This is for solid arrays, where there is no advection through 
    /// the array
    SolidArrayCV(SolidColumn),
}

impl From<SingleCVNode> for CVType {
    fn from(single_cv: SingleCVNode) -> Self {
        Self::SingleCV(single_cv)
    }
}

impl From<FluidArray> for CVType {
    fn from(fluid_array: FluidArray) -> Self {
        Self::FluidArrayCV(fluid_array)
    }
}

impl From<SolidColumn> for CVType {
    fn from(solid_array: SolidColumn) -> Self {
        Self::SolidArrayCV(solid_array)
    }
}

impl CVType {
    #[inline]
    /// gets the material 
    pub fn get_material(&mut self) -> Result<Material,ThermalHydraulicsLibError>{


        match self {
            CVType::SingleCV(single_cv_node) => {
                return Ok(single_cv_node.material_control_volume);
            },
            CVType::FluidArrayCV(fluid_array_cv) => {
                return Ok(fluid_array_cv.front_single_cv.material_control_volume);
            },
            CVType::SolidArrayCV(solid_array_cv) => {
                return Ok(solid_array_cv.front_single_cv.material_control_volume);
            },
        }
    }

    #[inline]
    fn get_temperature_vector(&mut self) -> 
    Result<Vec<ThermodynamicTemperature>,ThermalHydraulicsLibError>{
        match self {
            CVType::SingleCV(single_cv) => {
                let temperature = single_cv.get_temperature_from_enthalpy_and_set()?;

                let mut temp_vec: Vec<ThermodynamicTemperature> = vec![];

                temp_vec.push(temperature);

                return Ok(temp_vec);

            },
            CVType::FluidArrayCV(fluid_array_cv) => {
                return fluid_array_cv.get_temperature_vector();
            },
            CVType::SolidArrayCV(solid_array_cv) => {
                return solid_array_cv.get_temperature_vector();
            },

        }
    }
}
