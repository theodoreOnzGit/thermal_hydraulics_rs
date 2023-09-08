use crate::thermal_hydraulics_error::ThermalHydraulicsLibError;
use crate::heat_transfer_lib::control_volume_calculations::heat_transfer_entities::SingleCVNode;

use super::ArrayCVType;

impl ArrayCVType {
    /// returns a mutable borrow of the back cv
    #[inline]
    pub fn nested_back_cv_deref_mut(&mut self) -> 
    Result<&mut SingleCVNode, ThermalHydraulicsLibError> {

        let cv_ref = match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                &mut cartesian_array_cv.inner_single_cv
            },
            ArrayCVType::GenericPipe(fluid_arr) => {
                &mut fluid_arr.back_single_cv
            },

        };

        Ok(cv_ref)
    }
    /// returns a mutable borrow of the front cv
    #[inline]
    pub fn nested_front_cv_deref_mut(&mut self) -> 
    Result<&mut SingleCVNode, ThermalHydraulicsLibError> {

        let cv_ref = match self {
            ArrayCVType::Cartesian1D(cartesian_array_cv) => {
                &mut cartesian_array_cv.outer_single_cv
            },
            ArrayCVType::GenericPipe(fluid_arr) => {
                &mut fluid_arr.front_single_cv
            },
        };

        Ok(cv_ref)
    }
}
