# PWSDDE-cont
Periodic orbit continuation routines for piecewise-smooth delay differential equations
- automatic monodromy matrix formulation, stability analysis
- built in discontinuity induced bifurcation detection routines
- two parameter continuation of grazing and sliding bifurcations

## Installation
- Clone the repository or download the zip file containing the codes and place them in an arbitrary folder.
- Add this folder to the matlab path as: addpath(genpath("folder_name")).

## Usage
At first it is recommended to take a look at the three example continuation problems "fric_osc_demo.m", "bld_osc_demo.m", and "fric_osc_ODE_demo.m", which are included in the release. These can serve as great templates for creating PWSDDE-cont compatible continuation problems. 

- To define a continuation problem the users are expected to provide all necessary functions and Jacobians as seen in the "#demo" folder.
- Then a periodic orbit may be initialized using "sim_ns_dde.m" or manually, by providing the required fields.
- Periodic orbits can be corrected via Newton iteration using "orb_corr.m".
- Continuation of periodic orbits in on and two parameters are both done using "br12_cont.m".
- "orb_convert.m" may be used to change the solution signature of orbits found at bifurcation points.
- Finally a wide range of visualization options is available in the "plot tools" folder.

## References

The inner workings of the continuation codes and further information on the use of the code base is published in: [...]
