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
| `2D_cell_problem.py` | solver and visualization for 2D cell problem (mesh generation with `mshr`) | Figure 4 and $D^{\mathrm{eff}}$ in Section 5.1
| `3D_cell_problem.py` | solver and visualization for 3D cell problem (mesh generation with `mesh_gmsh.py`) |  Figure 5 and $D^{\mathrm{eff}}$ in Section 5.2
| `mesh_gmsh.py` | 3D mesh generation with Gmsh, called by the 3D script | ---
 

## Requirements
 
The code uses **legacy FEniCS (`dolfin`, `mshr`)**, which is no longer maintained.
 
| Package | Version | Needed for |
|---|---|---|
| Python | 3.9.16 | both |
| FEniCS (dolfin) | 2019.1.0 | both |
| mshr | 2019.1.0 | 2D script |
| Gmsh | 4.15.0 | 3D script |
| Mesa | 20.2.1 | 3D script (provides `libGLU`, required by the Gmsh Python package) |
| meshio | 5.3.4 | 3D script |
| h5py | 3.8.0 | 3D script |
| numpy | 1.24.3 | both |
| pandas | 2.0.1 | both |
| matplotlib | 3.7.1 | both |

## Usage
 
Set up the environment (tested on Chiméra/Sněhurka cluster: https://cluster.karlin.mff.cuni.cz/):
 
```bash
module load fenics/2019.1.0 mesa/20.2.1 meshio/5.3.4
```

The 3D script also needs the Gmsh Python package, which is not part of these modules. Install it with pip:
 
```bash
pip install --user gmsh==4.15.0
```
 
Run the 2D cell problem:
 
```bash
python 2D_cell_problem.py
```
 
Run the 3D cell problem:
 
```bash
python 3D_cell_problem.py
```
 
The 3D script generates its mesh itself (via `mesh_gmsh.py` and Gmsh), so no data needs to be downloaded.

### Settings used in the paper
 
Settings are constants at the top of each script. The values below reproduce the results in Sections 5.1 and 5.2 of the paper.
 
| Constant | Script | Meaning | Value |
|---|---|---|---|
| `MESH_RESOLUTION` | 2D | Mesh resolution passed to `mshr` | 100 |
| `FE_DEGREE` | 2D | Polynomial degree of the CG space | 3 |
| `MESH_TYPES` | 3D | Geometries to compute (0 to 3) | (0, 1, 2, 3)|
| `MESH_CHAR_LENGTH` | 3D | Maximal mesh size passed to Gmsh | 0.0125 |
| `FE_DEGREE` | 3D | Polynomial degree of the CG space | 1 |
| `SLICE_Z` | 3D | Height of the slice plot in Figure 5 | 0.5 |
 
The 2D script computes both of its geometries in one run (can be altered in `BLOOD_SINKS` setup). In the 3D script, select the geometries by editing `MESH_TYPES`, for example `MESH_TYPES = (0, 1)`.

## Output
 
Each run writes its results to one folder per geometry, inside `2Dmeshes/` (2D script) or `3Dmeshes/` (3D script).
 
2D:
 
```
2Dmeshes/
├── mesh0/
│   ├── mesh.png                                 mesh, neuron boundary blue, blood sink red
│   ├── solution_eta0.png, .xdmf, .h5            corrector for the first component
│   └── solution_eta1.png, .xdmf, .h5            corrector for the second component
└── mesh1/                                       same files for the second geometry
```
 
3D:
 
```
3Dmeshes/
├── mesh0/                                       one folder per geometry (mesh0 to mesh3)
│   ├── mesh0.png                                3D plot of the mesh boundary
│   ├── slice_z05_mesh0.png                      cross-section at z = 0.5
│   ├── mesh.xdmf, .h5                           mesh (ParaView)
│   ├── solution_eta0_on0.xdmf, .h5              corrector for the first component (ParaView)
│   ├── solution_eta0_scatter.png                quick-look plot of the corrector 
│   ├── solution_eta1_on0.*, solution_eta1_scatter.png    second component
│   ├── solution_eta2_on0.*, solution_eta2_scatter.png    third component
│   └── gmsh_mesh0.msh, gmsh_mesh0_mesh.*, gmsh_mesh0_mf.*, boundary_markers.*   files produced by Gmsh and meshio
├── mesh1/ ...
```
 
The `.xdmf` files can be opened in ParaView if the corresponding `.h5` is in the same folder.
 
The effective matrix is printed to the terminal at the end of each run, as a table and in LaTeX form.

 
## Citation
 
If you use this code, please cite the paper:
 
```bibtex
@article{BuJuRoWi,
  author  = {Bul{\'\i}{\v c}ek, M. and J{\"u}ngel, A. and Rognes, M. E. and Wintrov{\'a}, L.},
  title   = {Well-posedness and homogenization of a waste-protein clearance model in cortical tissue},
  doi = 
}
```

## License
MIT
