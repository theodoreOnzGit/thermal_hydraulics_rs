
# Patch Notes 

## v 0.0.13 

Moved the TUAS boussinesq solver to another 
[github repository](https://github.com/theodoreOnzGit/tuas_boussinesq_solver).
The thermal_hydraulics_rs solver now has the coupled natural circulation 
regression case as an example of how to use the tuas_boussinesq_solver in 
your own library. The prelude here also includes the tuas_boussinesq_solver.

## v 0.0.12 

Now, I tentatively call the boussinesq solver the tuas_boussinesq_solver,
which is named after the Tuas industrial area in Singapore.

It is also an acronym for Thermo-hydraulic Uniphase Advection and Convection 
Solvers (TUAS). 

In v0.0.12, I have also added all calibrated coupled DRACS loop results 
for datasets A, B and C within the SAM publication. It has matched the 
DRACS loop flowrate experimental data to within 6.1%, and pri loop flowrate 
experimental data to within 4.4%. In contrast to SAM, agreement with 
experimental data, the max error was 6.76% for the DRACS loop flowrate 
and 6.65% for the primary loop flowrate. See reference:

Zou, L., Hu, G., O'Grady, D., & Hu, R. (2021). Code validation of 
SAM using natural-circulation experimental data from the compact integral 
effects test (CIET) facility. Nuclear Engineering and Design, 377, 111144.

Given that the calibrated simulation with the TUAS solver agreed better with 
experimental data than SAM, I consider this validation effort successful.
The tests are parked under the pre_built_components module, where we have 
the ciet_steady_state_natural_circulation_test_components module. Inside that,
I put the coupled_dracs_loop_tests modules with dataset_a, dataset_b 
and dataset_c.


## v 0.0.11 

v0.0.11 has updated several dependencies including Peroxide.

Moreover, I have added uncalibrated coupled DRACS loop for CIET as a 
regression test case.

## v 0.0.10 

Added the Shell and Tube Heat Exchanger case by Du et al. and the isolated 
DRACS loop for CIET as test cases.

## v 0.0.9

Found out that compiling openblas on Linux Mint was not as straightforward 
as Arch Linux. Therefore I just opted to use the system libraries rather 
than build this from source. 

Also, updated all dependencies, tested on Arch Linux and Linux Mint.

## v 0.0.8

Sorted out bug for heater example as the wrong hydraulic diameter 
was used, so that the Nusselt correlation used was wrong. This caused 
the CIET heater temperature profile to be about 100 K hotter at all 
points across the whole heater. Changed the odd shaped pipe constructor 
to reduce occurrences of this bug. Also, changed 
the odd shaped pipe constructor for the FluidArray class so that 
hydraulic diameter must be given as an input. 

The main test now features a transient test with csv files printed.

Also, the example tests were updated. 

For neatness's sake, patch notes were moved into a separate file.

## v 0.0.7

Major speedup for given examples and SingleCVNode. This was profiled 
using flamegraph (somewhat easier to read than callgrind and 
valgrind) from the cargo repositories, licensed under MIT
or Apache License. 

SingleCVNode will 
now have a temperature attribute where temperature can be 
accessed many times per timestep, rather than repeatedly requiring 
calculation from specific enthalpy many times per timestep. This 
was severely bottlenecking down the simulation speed.

Another major bottleneck was the repeated construction of cubic splines 
when obtaining thermodynamic properties of steel whenever the 
function is accessed, especially thermal conductivity and 
heat capacity. After the function finishes returning the values,
the spline object is destroyed, and needs to be recreated. This 
was extremely computationally intensive. In general, repeated 
accessing of thermophysical properties multiple times during the 
timestep was bad for computation speed. It's crucial to calculate 
once and use the same stored value for each timestep.

## v 0.0.6 

Added some examples within source code, meant to show how 
to use the library to construct heat transfer components. 
These examples were tested to be relatively memory leak free. 
However, these are not quite optimised for best cpu power usage,
it could use some further work. Nevertheless, it serves as a good 
example, or tutorial on how to use the heat transfer library.
Along with the fluid component traits.
Copy and paste these to start your own. 

Also, corrected some bugs with connecting via advection 
HeatTransferInteractionType between array CVs.
