# Well-posedness and homogenization of a waste-protein clearance model in cortical tissue

This repository contains supplementary code for the paper
> Bulíček, M., Jüngel, A., Rognes, M. E. & Wintrová, L.
> Well-posedness and homogenization of a waste-protein clearance model in cortical tissue.

This code computes the effective diffusion matrix $D^{\mathrm{eff}}$ from the periodic cell problem for representative geometries in 2D and 3D. It reproduces results from Section 5 of the paper.


## Abstract
We formulate and analyze a mathematical model for the clearance of waste proteins from cortical tissue. The model tracks three species of waste proteins, distinguished by their size, and incorporates constant production at neuronal surfaces, reaction–diffusion in the bulk of interstitial fluid, and nonlinear dynamics on capillary surfaces. We establish the existence, boundedness, and uniqueness of a solution to the microscopic model. 
Using two-scale convergence, we derive a homogenized model and prove that the microscopic solutions converge to the corresponding macroscopic limit. The resulting model reveals how the microscopic arrangement and relative volume fractions of neurons and capillaries affect waste-protein transport and clearance at the tissue scale. Numerical simulations illustrate these effects and demonstrate how variations in cortical microstructure translate into differences in the macroscopic clearance dynamics.

## Contents
 
| File | Description | Result in the paper
|---|---|---|
| `2D_cell_problem.py` | 2D cell problem (mesh generation with `mshr`) | Figure 4 and $D^{\mathrm{eff}}$ in Section 5.1
| `3D_cell_problem.py` | 3D cell problem (solver, plots, output) |  Figure 5 and $D^{\mathrm{eff}}$ in Section 5.2
| `mesh_gmsh.py` | 3D mesh generation with Gmsh, called by the 3D script | ---
 

## Requirements
 
The code uses **legacy FEniCS (`dolfin`, `mshr`)**, which is no longer maintained.
 
| Package | Version | Needed for |
|---|---|---|
| Python | 3.9.16 | both |
| FEniCS (dolfin) | 2019.1.0 | both |
| mshr | 2019.1.0 | 2D script |
| Gmsh | 4.15.0 | 3D script |
| meshio | 5.3.4 | 3D script |
| h5py | 3.8.0 | 3D script |
| numpy | 1.24.3 | both |
| pandas | 2.0.1 | both |
| matplotlib | 3.7.1 | both |

## License
MIT
