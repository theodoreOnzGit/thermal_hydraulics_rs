

///
/// For this part, we determine how heat_transfer_entities interact 
/// with each other by using a function 
/// The function will take in a HeatTransferInteractionType enum
/// which you must first initiate
///
/// Notes for parallel computation: 
///
/// if we use mutable borrows early, it will be difficult to use 
/// this in a multithreaded manner
///
/// We can use interior mutability to allow for mutex locks and Arc 
/// types to be used so that parallel calculations can be performed 
///
/// lastly, we can just make sure we return a value from the placeholder 
/// function, and then manually map the heat transfer value to each 
/// heat_transfer_entities pair 
fn placeholder_function_suitable_for_parallel_computation(){}
