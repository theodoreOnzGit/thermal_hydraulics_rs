/// HeatTransferEntity module 
///
/// For practical reasons, using different functions to connect 
/// control volumes of various types (whether singleCV or arrayed control 
/// volumes) can be quite cumbersome 
///
/// To help the user connect these, I classify (and abstract) all control 
/// volumes and boundary conditions as HeatTransferEntity objects.
///
/// The basic use is that HeatTransferEntity objects are connected to 
/// each other by a user specified heat transfer interaction
/// 
pub mod heat_transfer_entities;

/// for fluid flow through pipes, these pipes will be represented by 
/// control volumes laterally coupled to one another. They come in 
/// insulated and non insulated forms
pub mod pipes;
