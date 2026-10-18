import os

import matplotlib.pyplot as plt
import mshr
import numpy as np
import pandas as pd
from dolfin import *

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
LX, LY = 1.25, 1.0                # cell dimensions
CELL_AREA = LX * LY               # full cell area, no cut-outs
MESH_RESOLUTION = 100             # resolution passed to mshr.generate_mesh
FE_DEGREE = 3                     # polynomial degree of the CG space

INNER_BOUNDARY = 2                # facet marker for neuron and sink walls
OUTER_BOUNDARY = 1                # facet marker for outer walls of the box (periodic)
NEURON_SINK_SPLIT_X = 0.98        # used for plotting and coloring the mesh; splits neurons (to the left) and sinks (to the right)

# Inclusions as ((x_min, y_min), (x_max, y_max))
NEURON = ((0.025, 0.025), (0.975, 0.975))
BLOOD_SINKS = {
    0: ((0.99375, 0.38125), (1.23125, 0.61875)),
    1: ((1.0875, 0.025), (1.146875, 0.975)),
}


# --------------------------------------------------------------------------
# Geometry and boundaries
# --------------------------------------------------------------------------
def generate_mesh(resolution, mesh_type):
    """Mesh of the cell with the neuron and the blood sink cut out."""
    print(f"Generating mesh type {mesh_type} with resolution {resolution}")
    geometry = mshr.Rectangle(Point(0.0, 0.0), Point(LX, LY))
    for (x0, y0), (x1, y1) in (NEURON, BLOOD_SINKS[mesh_type]):
        geometry = geometry - mshr.Rectangle(Point(x0, y0), Point(x1, y1))
    return mshr.generate_mesh(geometry, resolution)


class PeriodicBoundary(SubDomain):
    """Periodic identification of opposite edges of the cell."""
    def inside(self, x, on_boundary):
        # Master edges: left and bottom
        return bool((near(x[0], 0.0) or near(x[1], 0.0)) and on_boundary)

    def map(self, x, y):
        # Slave edges (right, top) are mapped onto the master edges
        y[0] = x[0] - LX if near(x[0], LX) else x[0]
        y[1] = x[1] - LY if near(x[1], LY) else x[1]


class InnerBoundary(SubDomain):
    """Identify inner walls: every boundary point not on the cell edge."""
    def inside(self, x, on_boundary):
        on_cell_edge = (near(x[0], 0.0) or near(x[0], LX)
                        or near(x[1], 0.0) or near(x[1], LY))
        return on_boundary and not on_cell_edge


# --------------------------------------------------------------------------
# Cell problem
# --------------------------------------------------------------------------
def solve_cell_problem(mesh, boundary_markers):
    """Return the effective matrix S and the list of correctors [eta_1, eta_2]."""
    V = FunctionSpace(mesh, "CG", FE_DEGREE, constrained_domain=PeriodicBoundary())
    dx = Measure("dx", domain=mesh)
    ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    normal = FacetNormal(mesh)

    u, v = TrialFunction(V), TestFunction(V)
    a = dot(grad(u), grad(v)) * dx

    S = np.zeros((2, 2))
    etas = []
    for i in range(2):
        L = -normal[i] * v * ds(INNER_BOUNDARY)   # Neumann data g = -n_i

        eta = Function(V)
        solve(a == L, eta)

        # eta is only defined up to a constant; fix it by zero mean of the DOF values. 
        # This affects the plots only, not S (which uses grad(eta)).
        eta.vector()[:] -= np.mean(eta.vector().get_local())
        etas.append(eta)

        for j in range(2):
            delta_ij = 1.0 if i == j else 0.0
            S[i, j] = assemble((delta_ij + grad(eta)[j]) * dx) / CELL_AREA
    return S, etas


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def output_folder(mesh_type):
    folder = f"2Dmeshes/mesh{mesh_type}"
    os.makedirs(folder, exist_ok=True)
    return folder


def save_mesh_plot(mesh, boundary_markers, mesh_type, filename="mesh.png"):
    """Plot the mesh as PNG; inner edges are blue (neuron) or red (blood sink)."""
    folder = output_folder(mesh_type)

    coords = mesh.coordinates()
    plt.figure(figsize=(8, 6))
    for facet in facets(mesh):
        edge = coords[facet.entities(0)]
        if boundary_markers[facet] != INNER_BOUNDARY:
            plt.plot(edge[:, 0], edge[:, 1], color="0.9", linewidth=1.0, zorder=1)
        else:
            color = "blue" if edge[:, 0].mean() < NEURON_SINK_SPLIT_X else "red"
            plt.plot(edge[:, 0], edge[:, 1], color=color, linewidth=1.0)
    plot(mesh, linewidth=0.3)
    plt.tight_layout()
    path = f"{folder}/{filename}"
    plt.savefig(path, dpi=300)
    plt.close()
    print(f"Mesh plot saved to {path}")


def save_corrector(eta, mesh_type, component):
    """Save a corrector as XDMF and as a PNG."""
    folder = output_folder(mesh_type)

    eta.rename(f"eta_{component}", f"eta_{component}")
    xdmf_path = f"{folder}/solution_eta{component}.xdmf"
    with XDMFFile(xdmf_path) as xf:
        xf.write(eta)

    plt.figure(figsize=(8, 6))
    plt.colorbar(plot(eta))
    plt.title(rf"Solution $\eta_{component}$ on Mesh {mesh_type}")
    plt.tight_layout()
    png_path = f"{folder}/solution_eta{component}.png"
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"Corrector eta{component} saved to {xdmf_path} and {png_path}")


def to_latex_matrix(A):
    lines = ["\\begin{bmatrix}"]
    lines += [" & ".join(f"{v:.3f}" for v in row) + " \\\\" for row in A]
    lines.append("\\end{bmatrix}")
    return "\n".join(lines)


def print_matrix(S, mesh_type):
    df = pd.DataFrame(S, columns=["∂x", "∂y"], index=["eta 1", "eta 2"])
    print(f"For mesh {mesh_type}, the effective diffusion matrix S is:")
    print(df.to_string(float_format="{:.3f}".format))
    print("In LaTeX:")
    print(to_latex_matrix(S))


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
for mesh_type in BLOOD_SINKS:
    mesh = generate_mesh(MESH_RESOLUTION, mesh_type)

    boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    InnerBoundary().mark(boundary_markers, INNER_BOUNDARY)
    PeriodicBoundary().mark(boundary_markers, OUTER_BOUNDARY)
    save_mesh_plot(mesh, boundary_markers, mesh_type)

    S, etas = solve_cell_problem(mesh, boundary_markers)
    for component, eta in enumerate(etas):
        save_corrector(eta, mesh_type, component)

    print_matrix(S, mesh_type)