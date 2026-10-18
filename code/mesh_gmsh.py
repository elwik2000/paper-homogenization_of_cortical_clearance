import gmsh
import meshio

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------
LX, LY, LZ = 1.25, 1.0, 1.0       # cell dimensions
TOL = 1e-6                        # tolerance for locating the outer faces

OUTER_TAG = 1                     # physical tag: outer (periodic) faces
INNER_TAG = 2                     # physical tag: inner (Neumann) faces
MATRIX_TAG = 1                    # physical tag: the volume
XDMF_DATA_NAME = "name_to_read"   # name of the marker data inside the XDMF files

# Outer faces as name: (axis, coordinate of the face)
OUTER_FACES = {
    "xmin": (0, 0.0), "xmax": (0, LX),
    "ymin": (1, 0.0), "ymax": (1, LY),
    "zmin": (2, 0.0), "zmax": (2, LZ),
}


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------
def _cut(occ, base, tools):
    """Subtract the volumes `tools` from the volume `base`; return the new volume tag."""
    result, _ = occ.cut([(3, base)], [(3, t) for t in tools],
                        removeObject=True, removeTool=True)
    return result[0][1]


def _concave_neuron(occ):
    """Cube [0.025, 0.975]^3 with a spherical dimple cut into each face."""
    x_0, size = 0.025, 0.95
    cube = occ.addBox(x_0, x_0, x_0, size, size, size)
    radius, depth = 0.7083, 0.1829
    offset = radius - depth              # sphere centre's distance from corresponding face
    lo, hi = x_0, x_0 + size        
    c = x_0 + size / 2                 # cube centre (each coordinate)
    centres = [(hi + offset, c, c), (lo - offset, c, c),
               (c, hi + offset, c), (c, lo - offset, c),
               (c, c, hi + offset), (c, c, lo - offset)]
    spheres = [occ.addSphere(x, y, z, radius) for x, y, z in centres]
    return _cut(occ, cube, spheres)


def _concave_sink(occ):
    """Square prism running along z, with a cylindrical dimple on each side."""
    x0, y0, size = 0.9925, 0.38, 0.24
    prism = occ.addBox(x0, y0, 0.0, size, size, LZ)
    radius, depth = 0.3585, 0.0207
    offset = radius - depth              # cylinder axis' distance from its side
    cx, cy = x0 + size / 2, y0 + size / 2
    axes = [(x0 + size + offset, cy), (x0 - offset, cy),
            (cx, y0 + size + offset), (cx, y0 - offset)]
    cylinders = [occ.addCylinder(x, y, 0.0, 0, 0, LZ, radius) for x, y in axes]
    return _cut(occ, prism, cylinders)


def _cylinder_sink(occ):
    return occ.addCylinder(1.1125, 0.5, 0.0, 0, 0, LZ, 0.11875)


def _build_inclusions(occ, mesh_type):
    """Create the neuron and the blood sink; return their volume tags."""
    if mesh_type == 0:                                                                      #flat geometry
        neuron = occ.addBox(0.11705, 0.11705, 0.11705, 0.7659, 0.7659, 0.7659)
        sink = occ.addBox(1.01, 0.395, 0.0, 0.2105, 0.2105, LZ)
    elif mesh_type == 1:                                                                    #convex geometry
        neuron = occ.addSphere(0.5, 0.5, 0.5, 0.475)
        sink = _cylinder_sink(occ)
    elif mesh_type == 2:                                                                    #concave geometry
        neuron = _concave_neuron(occ)
        sink = _concave_sink(occ)
    elif mesh_type == 3:                                                                    #mixed geometry
        neuron =  _concave_neuron(occ)
        sink = _cylinder_sink(occ)
    else:
        raise ValueError(f"unknown mesh_type {mesh_type}")
    return neuron, sink


def _classify_surfaces(domain_tag):
    """Split the boundary surfaces of the domain into the six outer faces
    (grouped by face) and the inner (neuron / sink) surfaces, using each
    surface's bounding box."""
    boundary = gmsh.model.getBoundary([(3, domain_tag)], oriented=False, combined=False)

    outer = {name: [] for name in OUTER_FACES}
    inner = []
    for dim, tag in boundary:
        bbox = gmsh.model.getBoundingBox(dim, tag)
        lower, upper = bbox[:3], bbox[3:]
        for name, (axis, position) in OUTER_FACES.items():
            is_flat = abs(upper[axis] - lower[axis]) < TOL
            if is_flat and abs(lower[axis] - position) < TOL:
                outer[name].append(tag)
                break
        else:
            inner.append(tag)
    return outer, inner


def _translation(dx, dy, dz):
    """4x4 affine translation matrix for periodicity, row-major, as expected by gmsh."""
    return [1, 0, 0, dx,  0, 1, 0, dy,  0, 0, 1, dz,  0, 0, 0, 1]


# --------------------------------------------------------------------------
# Mesh generation and conversion
# --------------------------------------------------------------------------
def generate_gmsh_mesh(mesh_type=0, char_length=0.05, out_prefix="mesh0"):
    """Generate the periodic tetrahedral mesh; return the names of the
    XDMF files holding the volume mesh and the facet markers."""
    gmsh.initialize()
    
    gmsh.model.add("cell")
    occ = gmsh.model.occ

    box = occ.addBox(0, 0, 0, LX, LY, LZ)
    neuron, sink = _build_inclusions(occ, mesh_type)
    domain_tag = _cut(occ, box, [neuron, sink])
    occ.synchronize()

    outer, inner = _classify_surfaces(domain_tag)
    for name, tags in outer.items():
        assert tags, f"no surface found on face {name} -- check geometry/tolerances"

    # Periodicity: slave face = master face + translation
    for slave, master, shift in (("xmax", "xmin", (LX, 0, 0)),
                                    ("ymax", "ymin", (0, LY, 0)),
                                    ("zmax", "zmin", (0, 0, LZ))):
        gmsh.model.mesh.setPeriodic(2, outer[slave], outer[master], _translation(*shift))

    # Physical groups (only these are written to the .msh file)
    gmsh.model.addPhysicalGroup(2, sum(outer.values(), []), tag=OUTER_TAG)
    gmsh.model.setPhysicalName(2, OUTER_TAG, "periodic")
    gmsh.model.addPhysicalGroup(2, inner, tag=INNER_TAG)
    gmsh.model.setPhysicalName(2, INNER_TAG, "neumann")
    gmsh.model.addPhysicalGroup(3, [domain_tag], tag=MATRIX_TAG)
    gmsh.model.setPhysicalName(3, MATRIX_TAG, "matrix")

    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", char_length)
    gmsh.model.mesh.generate(3)

    msh_path = f"{out_prefix}.msh"
    gmsh.write(msh_path)

    gmsh.finalize()

    _convert_to_xdmf(msh_path, out_prefix)
    return f"{out_prefix}_mesh.xdmf", f"{out_prefix}_mf.xdmf"


def _convert_to_xdmf(msh_path, out_prefix):
    """Split the .msh into a tetrahedra-only XDMF (volume mesh) and a
    triangles-only XDMF (facet markers), the format dolfin's XDMFFile expects."""
    msh = meshio.read(msh_path)

    tetra_mesh = meshio.Mesh(
        points=msh.points,
        cells={"tetra": msh.get_cells_type("tetra")},
        cell_data={XDMF_DATA_NAME: [msh.get_cell_data("gmsh:physical", "tetra")]},
    )
    triangle_mesh = meshio.Mesh(
        points=msh.points,
        cells={"triangle": msh.get_cells_type("triangle")},
        cell_data={XDMF_DATA_NAME: [msh.get_cell_data("gmsh:physical", "triangle")]},
    )

    meshio.write(f"{out_prefix}_mesh.xdmf", tetra_mesh)
    meshio.write(f"{out_prefix}_mf.xdmf", triangle_mesh)