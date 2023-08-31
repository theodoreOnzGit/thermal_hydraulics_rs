/// deals with fluid nodes in the core region
pub (in crate) mod core_fluid_node;

/// deals with fluid nodes as if they were in a shell region 
/// that means they are exposed to an inner region and an outer region
pub (in crate) mod shell_fluid_node;
