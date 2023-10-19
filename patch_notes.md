
# Patch Notes 

## v 0.0.9

Found out that compiling openblas on Linux Mint was not as straightforward 
as Arch Linux. Therefore I just opted to use the system libraries rather 
than build this from source. 

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
