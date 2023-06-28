## Courant Number Derivation for Heat Transfer

The courant number for heat transfer in 1D can be shown as follows:

$$Co = \frac{\alpha t_{timestep}}{L_{mesh}^2}$$

in 2D it can be [shown](https://skill-lync.com/student-projects/Analysis-of-Solution-Stability-in-a-2D-Heat-conduction-problem-08510)
as:


$$Co = \frac{\alpha t_{timestep}}{x_{mesh}^2} +
\frac{\alpha t_{timestep}}{y_{mesh}^2}$$

The maximum Courant number allowable here is not 1 but 0.25.


Now the courant number is a little harder to imagine for heat transfer,
unlike in fluid mechanics where we have:

$$Co_{fluid} = \frac{u \Delta t_{timestep}}{\Delta x}$$
or
$$Co_{fluid} = \frac{\dot{V} \Delta t_{timestep}}{V_{control\_volume}}$$

This is in the 1D case. 

Essentially, if the volumetric flowrate or mass flowrate is enough to
"flush" out all the control volume in one timestep, the Courant number
is too big.

For heat transfer, this intuition is not as obvious.

## Courant Number in releation to Fundamental Equations

### Conduction Boundary Conditions
Courant Number is used for transient heat transfer where in 1D:

$$m c_p \frac{\partial T}{\partial t} = 
-kA_{normal} \frac{\partial T}{\partial x}$$

This assumes that there is a control volume with one boundary of 
conduction heat transfer

$$\rho V_{control\_volume} c_p \frac{\partial T}{\partial t} = 
-kA_{normal} \frac{\partial T}{\partial x}$$

$$\frac{\partial T}{\partial t} = 
-\alpha \frac{A_{normal}}{V_{control\_volume}} \frac{\partial  T}{\partial x}$$



we note that the mesh size in x direction is:

$$x_{mesh} = \frac{V_{control\_volume}}{A_{normal}}$$

We could nondimensionalise the time with the timestep and x with
the mesh size in x direction

$$\frac{\partial T}{\partial t^*} = 
-\alpha \frac{t_{timestep}}{x_{mesh}^2} \frac{\partial  T}{\partial x^*}$$


The courant number then appears as:


$$Co = \frac{\alpha t_{timestep}}{x_{mesh}^2}$$

where:
$$x_{mesh} = \frac{V_{control\_volume}}{A_{normal}}$$

$$\frac{\partial T}{\partial t^*} = 
-Co \frac{\partial  T}{\partial x^*}$$

And it just so happens to have the same dimensions as the Fourier Number.

### Convection Boundary Condition
Now suppose we have a convection boundary condition, the Biot number comes
into play

$$m c_p \frac{\partial T}{\partial t} = 
-hA_{surface} (T_{surface} - T_{fluid})$$

$$\rho V_{control\_volume} c_p \frac{\partial T}{\partial t} = 
-hA_{surface} (T_{surface} - T_{fluid})$$


$$ \frac{\partial T}{\partial t} = 
\frac{-hA_{surface}}{\rho V_{control\_volume} c_p} (T_{surface} - T_{fluid})$$


we can define the lengthscale as volume to surface ratio:

$$L = \frac{V_{control\_volume}}{A_{surface}}$$


$$ \frac{\partial T}{\partial t} = 
\frac{-h}{\rho L c_p} (T_{surface} - T_{fluid})$$

Multiply top and bottom by conductivity of the control volume
($k_{solid}$),


$$ \frac{\partial T}{\partial t} = 
\frac{-h k}{\rho L k  c_p} (T_{surface} - T_{fluid})$$

$$ \frac{\partial T}{\partial t} = 
\frac{-h \alpha}{ L k  } (T_{surface} - T_{fluid})$$

nondimensionalising time and bringing out the fourier number,


$$ \frac{\partial T}{\partial t^*} = 
\frac{-h L \alpha t_{timescale}}{k L^2   } (T_{surface} - T_{fluid})$$

$$ \frac{\partial T}{\partial t^*} = 
-Bi Fo (T_{surface} - T_{fluid})$$

So in this case, the Courant number equivalent becomes the product
of
$$Co = Bi Fo$$

So how do we determine if something is table?

We can start with integration to get a sense of things

Let's integrate the dimensional form since it's more intuitive 
to integrate:
$$ \frac{d T}{d t} = 
\frac{-h \alpha}{ L k  } (T_{surface} - T_{fluid})$$

We assume $T_{surface} = T_{control\_volume}$ as for control volume,
lumped capacitance is always the case

$$ \frac{d (T_{cv} - T_{fluid})}{d t} = 
\frac{-h \alpha}{ L k  } (T_{cv} - T_{fluid})$$


$$ d\ \ln (T_{cv} - T_{fluid}) = 
\frac{-h \alpha}{ L k  } dt$$

we integrate from t =0  to t= $t_{elapsed}$


$$  (T_{cv}(t_{elapsed}) - T_{fluid})/(T_{cv\_initial} - T_{fluid})= 
\exp \left( \frac{-h \alpha}{ L k  } t_{elapsed} \right)$$

$$\theta = (T_{cv}(t_{elapsed}) - T_{fluid})/(T_{cv\_initial} - T_{fluid})$$
So the maximum temperature change that can happen is 

$$T_{cv\_initial} - T_{fluid}$$

So assuming that we have this change go to completion,


$$\exp \left(- \frac{h \alpha}{ L k  } t_{elapsed} \right) \approx 0$$

If we nondimensionalise 

$$t^* = \frac{t_{elapsed}}{t_{steady\_state}}$$

Where $t_{steady\_state}$ is the time taken for the heat transfer to reach 
completion, so that 

$$t^* = order\ of\ magnitude\ (1)$$

we can write,

$$ \exp \left(- Bi Fo\ t^* \right) \approx 0$$

I'll just define a threshold for steady state as:

$$\exp(- Bi\ Fo\ t^*) = 1 - 0.999 = 1e-3$$

take ln both sides,


$$- Bi\ Fo\ t^* = -6.9077$$

For a nice round easy to remember number,
$$ Bi\ Fo\ t^* \approx 7$$

This is the uppermost upper limit of the courant number defined
as (BiFo), as in the product of both cannot exceed 7, if not, the
heat transfer sort of "overshoots".

If the timestep is so big that it causes this, the simulation will surely
be unstable. If this Bi Fo becomes 7, just assume steady state happens instantly.

However, when it comes to timestepping, we know the gradient changes
at every timestep because the temperature difference becomes smaller.

So let's assume that the solver guesses the gradient of the change in 
$\theta$ with respect to time is linear, and for the timestep itself,
it assumes that the rate of change of temperature with respect to time
is constant. This would result in a linear relationship.

The gradient at t = 0 is simply -Bi Fo, 

the predicted evolution of temperature assuming a constant 
rate of change of temperature is:

$$\theta = 1 - Bi \ Fo \ t^*$$

If we were to plot the graphs of $\theta$ vs Bi Fo t*, and one graph
has exp(- Bi Fo t*) and the other has 1 - Bi Fo t*,

Then we would notice that their gradients diverge readily at Bi Fo = 0.25. 
So when Bi Fo  = 0.25, we would want the solver to re-evaluate the temperature
gradient to keep up with the changes. This is probably why we do not want the
courant number to be greater than 0.25. In fact Co = 0.1 to 0.2 is ideal so that
we can keep track of temperature changes accurately.

In the same way, Co should be 0.25 or less for the conduction heat flux case.
That's why Co should be 0.25 or less.

Again, to find Bi,

$$Bi = \frac{h_{fluid}L_{volume-to-surface-area-ratio}}{k_{solid}}$$

All other properties for $\alpha$ are solid thermophysical properties,
only the fluid heat transfer coefficient is a property of the fluid.

### Transport Boundary Conditions

Now suppose that we have fluid flowing in and out of the control volume
bringing in and out the scalar, in this case, temperature or more correctly,
enthalpy

The [scalar transport equation](https://www.cfd-online.com/Wiki/Generic_scalar_transport_equation) can be written as follows
for enthalpy:

$$  \frac{\partial \rho h}{\partial t} + \nabla \bullet (\rho \vec{u} h) =
\nabla \bullet (\Gamma \nabla h )$$

$\Gamma$ is the diffusion coefficient for the scalar, in this case it is enthalpy.
It will be somewhat different from the diffusion coefficient for temperature.

Now, if we just consider one control volume:

$$m_{cv} c_p \frac{\partial T}{\partial t} = \dot{m}_{in}h_{in} -
\dot{m}_{out} h_{cv}$$

we can simplify the terms assuming pressures don't change so much that
the $c_p$ changes:

$$dh = c_p dT$$

$$m_{cv} \frac{\partial h_{cv}}{\partial t} = \dot{m}_{in}h_{in} -
\dot{m}_{out} h_{cv}$$

For this case, the Courant Number is quite straightforward (assuming conduction
from any of those prior control volumes is quite negligible).

The momentum transport 
Courant number here will then be the OpenFOAM expression

$$Co = 0.5 \frac{\dot{m}_{in} t_{timestep}}{m_{cv}}
+ 0.5 \frac{\dot{m}_{out} t_{timestep}}{m_{cv}}$$

However, if we were to once again perform our analysis using

$$dh_{cv} = d(h_{cv} - h_{in})$$

and for a simple pipe:

$$\dot{m}_{in} = \dot{m}_{out}$$


$$m_{cv} \frac{\partial (h_{cv} - h_{in})}{\partial t} = -\dot{m}_{in}(h_{cv} - 
h_{in})$$

we get a similar expression for the courant number...

$$m_{cv} \frac{\partial \ln (h_{cv} - h_{in})}{\partial t} = 
-\dot{m}_{in}$$

$$\ln(h_{cv} - h_{in})_{t=t} - \ln(h_{cv}-h_{in})_{t=0} 
= - \frac{\dot{m}_{in}}{m_{cv}} t$$

$$(h_{cv} - h_{in}) /(h_{cv}-h_{in})_{t=0} 
= \exp \left(- \frac{\dot{m}_{in}}{m_{cv}} t \right)$$

We see now the courant number takes the form:

$$Co = \frac{\dot{m}_{in} t_{timestep}}{m_{cv}}$$

and in general, for a control volume with flows in and out,

$$Co = 0.5 \sum_i \frac{|\dot{m}_{i}| t_{timestep}}{m_{cv}} $$

This is taken from OpenFOAM's expression of Courant Number.

In this case though, since we are concerned about scalar transport,
the cournat number limit for this form is not 1, but instead 0.25.


### What if Courant Number Exceeds 0.25?

It is evident that the Courant number for heat transfer will be the main
bottleneck in timestep for most high prandtl number fluids as the control
volume Courant number for heat transfer due to flows is already
more stringent a requirement than the courant number for fluid flow because
both use the same expression.

$$Co = 0.5 \sum_i \frac{|\dot{m}_{i}| t_{timestep}}{m_{cv}} $$

So what if we exceed the courant number? Basically, we cut the 
timestep down locally. Hopefully rust is fast enough for this not to be an
issue.

Here is some rust code that could depict what could happen if courant 
number exceeds 0.25

```rust
// .... (this code has not been tested)

let Co: f64 = Co_inflow_and_outflow + Co_conduction + Co_convection

// declare a variable first

let mut timestep_intervals: u16;
let mut timestep_to_use: f64;

if Co > 0.25 {
    // check how many times the Co exceeds 0.25

    let times_exceeded: f64 = Co/0.25;

    // probably need to convert here to u16
    // in actual rust code
    timestep_intervals = times_exceeded.ceil();

    timestep_to_use = user_supplied_timestep/timestep_intervals;
} else {
    // if Co < 0.25 then timestep_intervals = 1
    timestep_intervals = 1;
    timestep_to_use = user_supplied_timestep;

}

for i in 0..timestep_intervals {
    // run calculations here

}

// some cleaning up ..

return control_volume_temperature;


```


## misc thoughts
### Conduction in the Fluid across adjacent control volumes

In this analysis of course, we assumed that the axial conduction in the 
direction of flow was negligible compared to the enthalpy transport due
to forced flow.

This is quite true except when flows are extremely slow moving, or we have
some low Prandtl Number fluid.

For high Prandtl number fluids such as oil, the momentum diffusivity is several
times more than thermal diffusivity. Therefore advection tends to happen
many times faster than thermal conduction.

However, this conduction is important in the case of blocked flow, or flows
where the temperature gradient is extremely high such as near the heater or
heat exchangers. 

In general for pipes, temperature gradients are small enough compared to
mass flowrates that one can neglect the conduction effects compared to
advection.













