/// For CIET natural circulation tests,
///
/// there are two sets of tests.
///
/// The first one isolates the DRACS 
///
/// in this test, the DRACS loop with zero parasitic heat loss is 
/// compared against the analytical solution by Scarlat for a 
/// natural circulation loop with zero heat loss
///
/// Three sets of tests are performed. 
/// One with TCHX outlet temperature of 46 C 
/// Second with TCHX outlet temperature of 40 C
/// Third with TCHX outlet temperature of 35 C
///
/// These are enforced as part of the boundary conditions of the TCHX
///
/// For this set of tests, we do not worry about real-time calculations 
/// just yet
///
/// 
pub mod zero_parasitic_heat_loss_isolated_dracs_loop_tests;
