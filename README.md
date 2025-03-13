# PWSDDE-cont
Periodic orbit continuation routines for piecewise-smooth delay differential equations (PWS-DDEs)
- Chebyshev collocation based numeric solution of the governing multi-point boundary value problem (MP-BVP)
- automatic monodromy matrix formulation, stability analysis, and detection of stability related bifurcations
- built in discontinuity induced bifurcation detection routines (grazing, sliding)
- two parameter continuation of grazing, sliding, or user defined bifurcations
- [COCO](https://sourceforge.net/projects/cocotools/) *"core"* compatible definition of the governing MP-BVP and its monitor functions
- limited treatment of neutral delays (WARNING: forward propagation of discontinuities may lead to major interpolation errors!)
- partial support for state dependent delays (only explicit dependence on the present state is allowed!)

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=Tudesz/PWSDDE-cont&file=bld_osc_demo.m)

## Installation
- Clone the repository or download the zip file containing the codes and place them in an arbitrary folder.
- Add this folder to the matlab path as: addpath(genpath("<folder_name>")).
- To make use of COCO when running the included demo codes, replace "<COCO_dir>" with the appropriate path.

## Usage
At first it is recommended to take a look at the four example continuation problems: **bld_osc_demo.m**, **fric_osc_demo.m**, **nlne_rob_demo.m**, and **impd_osc_demo.m**, which are included in the release. These can serve as great templates for creating PWSDDE-cont compatible continuation problems. A slightly different formalism for problems with state dependent delays is illustrated by **sd_broach_demo.m**. The most important routines of the codebase are all found in the *pwsdde cont* folder. For most continuation problems the use of these algorithms should be sufficient. For more intricate tasks, however, it is also worthwhile to take a deeper look into the *_toolbox* folder.

1) To define a continuation problem, the users are expected to provide all necessary functions and Jacobians as demonstrated by the examples found in the *_system def* folder. These should be coded up in three (four) separate functions:
    - **f(x,xd,p,mode,type,l)**: right hand side of the PWS-DDE and its Jacobians
    - **e(x,xd,p,id,type,l)**: event functions, event maps, mode transitions, and corresponding Jacobians
    - **tau(p,ind,type)**: time delays and their parameter Jacobians
    - **p(t,x,xd,p)**: an optional function for better visualization of the found periodic orbits via **plot_orb.m**
2) A periodic orbit may be initialized automatically from transient simulations using **sim_ns_dde.m**, or manually, by providing the required fields:
    - *sig*: solution signature stored as an event list (contains *id* identifiers for the **e( )** function)
    - *M*: Chebyshev mesh size (recommended 20)
    - *n*: state dimension
    - *p*: parameter vector
    - *U*: state variable vector (discrete solution on the piecewise-Chebyshev mesh)
    - *T*: segment lengths (times between events)
3) Periodic orbits can be corrected through the solution of the governing MP-BVP via **orb_corr.m**.
    - passing an additional *bifs* field also allows the identification of grazing, sliding, or user defined bifurcation points
    - the stability of the corrected orbits may be assesed via **orb_stab.m**
    - in cases of changes in the solution signature (*sig*), orbits may be projected to a new piecewise-mesh using **orb_convert.m**
4)  Periodic orbits can be followed either via **br12_cont_fix.m** or **br12_cont_adapt.m** employing a fixed-step or an adaptive pseudo-arclength method.
    - before calling either function an *opts* structure should be initialized using **br12_opts.m**, which contains *pi*, the index of the active continuation parameter in *p* 
    - for two parameter continuation *pi* should include two indices, and passing the *bifs* field is mandatory
    - the branch data structure returned by these functions contains information on encountered special points under the *bif_type* field
    - for further continuation options, the users are referred to the headers of **br12_cont_fix.m** and **br12_cont_adapt.m**
5) Finally, a wide range of visualization options is available in the *plot tools* folder.
    - **plot_res.m** may be used to visualize transient simulation results obtained via **sim_nds_dde.m**
    - **plot_orb.m** is used for visualizing the found periodic orbits (can be customized by defining a **p( )** function)
    - **plot_spectrum.m** can visualize Floquet multipliers on the Re-Im plane
    - **plot_br1_ampl** displays the amplitude of periodic orbits in a branch with respect to one bifurcation parameter
    - **plot_br1_norm** plots the norm of the MP-BVP solution for all elements of a branch with respect to one bifurcation parameter
    - **plot_br2_par** plots one bifurcation parameter with respect to another one
    - **plot_br2_3D** can be used to visualize either the norm or amplitude of all orbits in a branch with respect to two bifurcation parameters
    - **anim_br_orb.m** creates an animation by executing **plot_orb.m** for all orbits in a branch
    - **anim_br_spectrum.m** animates the progression of Floquet multipliers along a branch of periodic orbits

## References

The inner workings of the continuation algorithms and further information on the use of the codebase is published in the article
[Bifurcation analysis of piecewise-smooth engineering systems with delays through numeric continuation of periodic orbits](https://doi.org/10.1007/s11071-024-10188-8). If the present framework is used for research purposes, please include a reference to this work.

