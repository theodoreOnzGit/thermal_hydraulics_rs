/// Example 1: 
///
/// This example shows how to create a simple pipe
/// using the FluidComponent and FluidPipeCalcPressureLoss,
/// traits
///
/// this is by no means the best way to do it, but its a start
/// remember to use the relevant imports in the fluid component
/// tests
///
/// it is made of copper, 1m long, 2 in in diameter
///
/// This does not take inclined angles into consideration yet
pub mod air_pipe_example;


/// Example 2:
///
/// We saw previously how to create an air pipe
/// now we shall make a slanted water pipe
/// with some internal pressure source (as if it had a pump attached
/// to it)
///
/// we shall improve on how we can create the pipes
/// to do so, we shall use the FluidComponent trait and the 
/// FluidPipeCalcPressureChange trait
///
pub mod water_pipe_example;