import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from dolfin import *
from mesh_gmsh import (LX, LY, LZ, INNER_TAG, OUTER_TAG, XDMF_DATA_NAME,
                       generate_gmsh_mesh)

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
MESH_TYPES = (0, 1, 2, 3)         # geometries to run (0-3, see mesh_gmsh.py)
MESH_CHAR_LENGTH = 0.0125          # maximal mesh size passed to gmsh
FE_DEGREE = 1                     # polynomial degree of the CG space

CELL_VOLUME = LX * LY * LZ        # full cell volume, no cut-outs
SLICE_Z = 0.5                     # used for plotting the mesh; height of the slice plot
NEURON_SINK_SPLIT_X = 0.985       # used for plotting and coloring the mesh; splits neurons (to the left) and sinks (to the right)


# --------------------------------------------------------------------------
# Periodic boundary
# --------------------------------------------------------------------------
class PeriodicBoundary(SubDomain):
    """Periodic identification of opposite faces of the cell."""

    def inside(self, x, on_boundary):
        # Master faces: x = 0, y = 0, z = 0
        return bool((near(x[0], 0.0) or near(x[1], 0.0) or near(x[2], 0.0)) and on_boundary)

    def map(self, x, y):
        # Slave faces (x = LX, y = LY, z = LZ) are shifted onto the master faces,
        y[0] = x[0] - LX if near(x[0], LX) else x[0]
        y[1] = x[1] - LY if near(x[1], LY) else x[1]
        y[2] = x[2] - LZ if near(x[2], LZ) else x[2]


# --------------------------------------------------------------------------
# Mesh loading and checks
# --------------------------------------------------------------------------
def load_mesh(mesh_xdmf, markers_xdmf):
    """Read the volume mesh and the facet markers written by mesh_gmsh.py."""
    mesh = Mesh()
    with XDMFFile(mesh_xdmf) as infile:
        infile.read(mesh)

    mvc = MeshValueCollection("size_t", mesh, mesh.topology().dim() - 1)
    with XDMFFile(markers_xdmf) as infile:
        infile.read(mvc, XDMF_DATA_NAME)
    return mesh, cpp.mesh.MeshFunctionSizet(mesh, mvc)


def check_periodicity(mesh, tol=1e-6):
    """Sanity check: every vertex on a master face (x, y or z = 0) has a
    slave vertex on the opposite face."""
    coords = mesh.coordinates()
    for axis, length in enumerate((LX, LY, LZ)):
        master = coords[np.isclose(coords[:, axis], 0.0, atol=tol)]
        slave = coords[np.isclose(coords[:, axis], length, atol=tol)]
        unmatched = 0
        for vertex in master:
            target = vertex.copy()
            target[axis] = length
            if np.linalg.norm(slave - target, axis=1).min() > tol:
                unmatched += 1
        print(f"Periodicity in {'xyz'[axis]}: {unmatched} of {len(master)} master vertices unmatched ({len(slave)} on the opposite face)")


# --------------------------------------------------------------------------
# Cell problem
# --------------------------------------------------------------------------
def solve_cell_problem(mesh, boundary_markers):
    """Return the effective matrix S and the list of correctors [eta_1, eta_2, eta_3]."""
    V = FunctionSpace(mesh, "CG", FE_DEGREE, constrained_domain=PeriodicBoundary())
    print(f"Function space: {V.dim()} degrees of freedom")
    dx = Measure("dx", domain=mesh)
    ds = Measure("ds", domain=mesh, subdomain_data=boundary_markers)
    normal = FacetNormal(mesh)

    u, v = TrialFunction(V), TestFunction(V)
    a = dot(grad(u), grad(v)) * dx

    S = np.zeros((3, 3))
    etas = []
    for i in range(3):
        L = -normal[i] * v * ds(INNER_TAG)        # Neumann data g = -n_i

        eta = Function(V)
        solve(a == L, eta, solver_parameters={"linear_solver": "mumps"})

        # eta is only defined up to a constant; fix it by zero mean of the DOF values. 
        # This affects the plots only, not S (which uses grad(eta)).
        eta.vector()[:] -= np.mean(eta.vector().get_local())
        etas.append(eta)

        for j in range(3):
            delta_ij = 1.0 if i == j else 0.0
            S[i, j] = assemble((delta_ij + grad(eta)[j]) * dx) / CELL_VOLUME
    return S, etas


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------
def output_folder(mesh_type):
    folder = f"3Dmeshes/mesh{mesh_type}"
    os.makedirs(folder, exist_ok=True)
    return folder


def inner_facet_color(triangle):
    """Blue (neuron) or red (blood sink), judged by the facet's centroid."""
    return "blue" if triangle[:, 0].mean() < NEURON_SINK_SPLIT_X else "red"


def save_mesh_plot(mesh, boundary_markers, mesh_type, filename):
    """3D plot of the boundary facets: outer faces white, neuron blue, sink red."""
    folder = output_folder(mesh_type)
    coords = mesh.coordinates()
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot(111, projection="3d")

    for facet in facets(mesh):
        marker = boundary_markers[facet]
        triangle = coords[facet.entities(0)]      # the 3 vertices of the facet
        if marker == OUTER_TAG:
            color = "white"
        elif marker == INNER_TAG:
            color = inner_facet_color(triangle)
        else:
            continue
        poly = Poly3DCollection([triangle], alpha=0.5)
        poly.set_facecolor(color)
        poly.set_edgecolor(color)
        ax.add_collection3d(poly)

    ax.set_xlim(coords[:, 0].min(), coords[:, 0].max())
    ax.set_ylim(coords[:, 1].min(), coords[:, 1].max())
    ax.set_zlim(coords[:, 2].min(), coords[:, 2].max())
    ax.set_box_aspect((LX, LY, LZ))

    plt.tight_layout()
    path = f"{folder}/{filename}"
    plt.savefig(path, dpi=300)
    plt.close()
    print(f"Mesh plot saved to {path}")


def save_mesh_xdmf(mesh, boundary_markers, mesh_type):
    """Save mesh and facet markers for ParaView."""
    folder = output_folder(mesh_type)
    with XDMFFile(f"{folder}/mesh.xdmf") as xf:
        xf.write(mesh)
    with XDMFFile(f"{folder}/boundary_markers.xdmf") as xf:
        xf.write(boundary_markers)
    print(f"Mesh and facet markers saved to {folder}")


def triangle_plane_intersection(triangle, z_plane):
    """End points of the intersection of a triangle with the plane z = z_plane."""
    dist = triangle[:, 2] - z_plane
    side = np.where(np.abs(dist) < 1e-6, 0, np.sign(dist))

    points = []
    for k in range(3):
        m = (k + 1) % 3
        if side[k] == 0:
            points.append(triangle[k])                          # vertex on the plane
        elif side[k] * side[m] < 0:                             # edge crosses the plane
            t = dist[k] / (dist[k] - dist[m])
            points.append(triangle[k] + t * (triangle[m] - triangle[k]))
    return points if len(points) == 2 else None


def save_slice_plot(mesh, boundary_markers, mesh_type, z_slice, filename):
    """2D plot of the inner boundaries cut by the plane z = z_slice."""
    coords = mesh.coordinates()
    fig, ax = plt.subplots(figsize=(8, 6))

    for facet in facets(mesh):
        if boundary_markers[facet] != INNER_TAG:
            continue
        triangle = coords[facet.entities(0)]
        segment = triangle_plane_intersection(triangle, z_slice)
        if segment is not None:
            segment = np.array(segment)
            ax.plot(segment[:, 0], segment[:, 1],
                    color=inner_facet_color(triangle), linewidth=1.5)

    ax.set_xlim(coords[:, 0].min(), coords[:, 0].max())
    ax.set_ylim(coords[:, 1].min(), coords[:, 1].max())
    ax.set_aspect("equal")

    plt.tight_layout()
    path = f"{output_folder(mesh_type)}/{filename}"
    plt.savefig(path, dpi=300)
    plt.close()
    print(f"Slice plot saved to {path}")


def save_corrector(mesh, eta, mesh_type, component):
    """Save a corrector as XDMF and as a sanity-check scatter plot."""
    folder = output_folder(mesh_type)

    eta.rename(f"eta_{component}", f"eta_{component}")
    xdmf_path = f"{folder}/solution_eta{component}_on{mesh_type}.xdmf"
    with XDMFFile(xdmf_path) as xf:
        xf.write(eta)

    # Vertex values colored by value: a sanity check without leaving Python
    coords = mesh.coordinates()
    values = eta.compute_vertex_values(mesh)
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=values, cmap="viridis", s=6)
    fig.colorbar(sc, ax=ax, shrink=0.7)
    ax.set_title(rf"Solution $\eta_{component}$ on Mesh {mesh_type} (vertex values)")
    plt.tight_layout()
    png_path = f"{folder}/solution_eta{component}_scatter.png"
    plt.savefig(png_path, dpi=300)
    plt.close()
    print(f"Corrector saved to {xdmf_path} and {png_path}")


def to_latex_matrix(A):
    lines = ["\\begin{bmatrix}"]
    lines += [" & ".join(f"{v:.3f}" for v in row) + " \\\\" for row in A]
    lines.append("\\end{bmatrix}")
    return "\n".join(lines)


def print_matrix(S, mesh_type):
    df = pd.DataFrame(S, columns=["∂x", "∂y", "∂z"], index=["eta 1", "eta 2", "eta 3"])
    print(f"For mesh {mesh_type}, the effective diffusion matrix S is:")
    print(df.to_string(float_format="{:.3f}".format))
    print("In LaTeX:")
    print(to_latex_matrix(S))


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------
for mesh_type in MESH_TYPES:
    mesh_xdmf, markers_xdmf = generate_gmsh_mesh(
        mesh_type=mesh_type, char_length=MESH_CHAR_LENGTH,
        out_prefix=f"{output_folder(mesh_type)}/gmsh_mesh{mesh_type}")
    mesh, boundary_markers = load_mesh(mesh_xdmf, markers_xdmf)
    print(f"Mesh {mesh_type}: {mesh.num_cells()} cells, hmin = {mesh.hmin():.4g}, hmax = {mesh.hmax():.4g}")

    save_mesh_plot(mesh, boundary_markers, mesh_type, f"mesh{mesh_type}.png")
    save_mesh_xdmf(mesh, boundary_markers, mesh_type)
    save_slice_plot(mesh, boundary_markers, mesh_type, SLICE_Z, f"slice_z05_mesh{mesh_type}.png")
    check_periodicity(mesh)

    S, etas = solve_cell_problem(mesh, boundary_markers)
    for component, eta in enumerate(etas):
        save_corrector(mesh, eta, mesh_type, component)

    print_matrix(S, mesh_type)
