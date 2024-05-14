
use super::one_d_fluid_array_with_lateral_coupling::FluidArray;


/// this is essentially a parallel collection of 
/// identical 1D FluidArray structs 
/// 
/// This is meant for heat exchangers which contain parallel collections 
/// of tubes which are represented by one single pipe for the 
/// sake of calculation convenience
///
#[derive(Debug,Clone,PartialEq)]
pub struct ParallelFluidArray {

    /// this is the representative one dimensional pipe/fluid_component
    /// which we use to model one of the many 
    /// parallel pipes/fluid_components
    /// within this parallel fluid array
    /// 
    pub representative_fluid_component: FluidArray,

    /// number of parallel pipes/fluid_components 
    pub number_of_parallel_tubes: u32,



}
