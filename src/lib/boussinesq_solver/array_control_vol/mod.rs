
/// contains matrix calculations specific to fluid nodes
/// arranged in a 1D array
///
/// These are standalone, and not abstracted under an arrayCV
/// struct, they only use SingleCVNode structs and representative 
/// arrays to represent an array of control volumes
pub mod standalone_fluid_nodes;

/// contains matrix calculations specific to solid nodes 
/// these are meant to represent the "shell" of the pipe 
/// or any kind of solid material in the pipe
///
/// These are standalone, and not abstracted under an arrayCV
/// struct, they only use SingleCVNode structs and representative 
/// arrays to represent an array of control volumes
pub mod standalone_solid_nodes;


pub mod one_dimension_cartesian_conducting_medium;
