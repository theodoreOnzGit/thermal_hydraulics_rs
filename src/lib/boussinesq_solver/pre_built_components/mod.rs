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

/// for fluid flow through non insulated pipes and fluid components, these pipes will 
/// be represented by control volumes laterally coupled to one 
/// another. 
pub mod non_insulated_fluid_components;


/// for fluid flow through insulated and fluid components with one layer 
/// of insulation, these pipes will 
/// be represented by control volumes laterally coupled to one 
/// another. 
pub mod insulated_pipes_and_fluid_components;

/// represents one dimensional solid structure
pub mod one_d_solid_structure;

/// represents the old CIET heater version 2 based on 
/// https://escholarship.org/uc/item/0362h3zf
///
/// Ong, T. K. C. (2024). Digital Twins as 
/// Testbeds for Iterative Simulated Neutronics Feedback 
/// Controller Development (Doctoral dissertation, UC Berkeley).
pub mod ciet_heater_version_2_bare;

/// represents the old CIET struct support codes based on
/// https://escholarship.org/uc/item/0362h3zf
///
/// Ong, T. K. C. (2024). Digital Twins as 
/// Testbeds for Iterative Simulated Neutronics Feedback 
/// Controller Development (Doctoral dissertation, UC Berkeley).
pub mod ciet_struct_supports;

/// represents the CIET heater top and bottom head codes based on
/// https://escholarship.org/uc/item/0362h3zf
///
/// Ong, T. K. C. (2024). Digital Twins as 
/// Testbeds for Iterative Simulated Neutronics Feedback 
/// Controller Development (Doctoral dissertation, UC Berkeley).
pub mod ciet_heater_top_and_bottom_head_bare;


/// represents the CIET static mixer 10 codes based on
/// https://escholarship.org/uc/item/0362h3zf
///
/// Ong, T. K. C. (2024). Digital Twins as 
/// Testbeds for Iterative Simulated Neutronics Feedback 
/// Controller Development (Doctoral dissertation, UC Berkeley).
pub mod ciet_static_mixer_mx_10;

/// ciet components for pipes and valves for use in the isothermal test
///
/// Zweibaum, N. (2015). Experimental validation of passive safety 
/// system models: Application to design and optimization of 
/// fluoride-salt-cooled, high-temperature reactors. University 
/// of California, Berkeley.
///
/// In my master's thesis, heat structure information was not included. However, 
/// I shall include them in this round
pub mod ciet_isothermal_test_components;
