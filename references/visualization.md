# Quantum Chemistry Visualization in Python

A practical guide covering orbital visualization, potential energy surfaces, reaction paths, and correlation diagrams.

---

## 1. Orbital Visualization with orbital-viz

[orbital-viz](https://pypi.org/project/orbital-viz/) ([GitHub](https://github.com/AleksandraLeszczyk/OrbitalViz)) is a Python library for interactive visualization of molecular orbitals. It evaluates contracted GTO basis functions on a 3D grid and renders ψ isosurfaces as interactive Plotly figures. Released June 2026, requires Python ≥ 3.11.

**Install:**
```bash
pip install orbital_viz
```

**Dependencies:** NumPy, SciPy, Plotly (no heavy QC code required at runtime).

### From a Molden file

The fastest path: parse geometry, basis, and MO coefficients directly from a standard Molden file.

```python
from orbital_viz import parse_molden_to_dict, BasisGTO, read_molden_c_matrix, plot_molecular_orbital

# Parse geometry and basis set
data = parse_molden_to_dict("molecule.molden")
basis = BasisGTO(**data)

# Read MO coefficient matrix  →  shape (n_basis, n_mo)
C = read_molden_c_matrix("molecule.molden", basis.n_basis)

# Plot MO index 5 (zero-based) – returns a Plotly Figure
fig = plot_molecular_orbital(basis, C, n=5)
fig.show()                         # interactive in Jupyter / browser
fig.write_html("mo5.html")         # save as standalone HTML
```

### From a custom basis (e.g. PyBEST output)

Build `BasisGTO` manually when your coefficients come from a non-Molden source.

```python
from orbital_viz import BasisGTO, plot_molecular_orbital
import numpy as np

basis = BasisGTO(
    atoms=[6, 6],                                    # atomic numbers (C, C)
    coordinates=[[-0.67, 0, 0], [0.67, 0, 0]],      # positions in Å
    number_of_primitives=[3, 3],
    contraction=[0.154, 0.535, 0.444, 0.154, 0.535, 0.444],
    Alpha=[71.6, 13.0, 3.53, 71.6, 13.0, 3.53],
    shell_types=[0, 0],          # 0 = s-type shell
    shell_to_atom=[0, 1],
)

C = np.load("mo_coefficients.npy")    # shape (n_basis, n_mo)

# HOMO on dark background – great for presentation figures
fig = plot_molecular_orbital(basis, C, n=4, isovalue=0.04, dark_bg=True)
fig.show()
```

### API reference

**`parse_molden_to_dict(filepath)`** — Parses a Molden file and returns a dict with keys `atoms`, `coordinates`, `number_of_primitives`, `contraction`, `Alpha`, `shell_types`, `shell_to_atom`. Pass directly to `BasisGTO(**data)`.

**`read_molden_c_matrix(filepath, n_basis)`** — Returns `np.ndarray` of shape `(n_basis, n_mo)`.

**`BasisGTO`** — Contracted GTO basis set. Handles Cartesian and spherical (solid-harmonic) AOs up to arbitrary angular momentum. Key attributes: `n_atoms`, `n_basis`.

**`plot_molecular_orbital(basis, mo_coeffs, n, **kwargs)`** — Core rendering function. Returns a `plotly.graph_objects.Figure`.

| Parameter | Default | Description |
|---|---|---|
| `grid_points` | `50` | Grid resolution per axis (increase for smoother surfaces) |
| `isovalue` | `0.05` | ψ value at which isosurfaces are drawn |
| `padding` | `2.0` | Extra space (Å) around the molecular bounding box |
| `opacity` | `0.45` | Isosurface opacity (0–1) |
| `atom_scale` | `0.25` | Scale factor on VDW radii for atom spheres |
| `bond_radius` | `0.08` | Bond cylinder radius (Å) |
| `title` | `None` | Figure title; defaults to `"MO #n"` |
| `dark_bg` | `False` | Dark background (good for glowing orbital renders) |
| `show_labels` | `True` | Show atom labels |

Example notebooks are in the [examples/](https://github.com/AleksandraLeszczyk/OrbitalViz/tree/main/examples) directory: `ethylene_from_pybest.ipynb` and `cyclobutadiene_from_molden.ipynb`.

### Other orbital visualization tools

For wave function post-processing from Gaussian/ORCA/MOLPRO output (cube files, density grids), [ORBKIT](https://orbkit.github.io/) remains a good complement to orbital-viz. For periodic systems, [orbvis](https://pypi.org/project/orbvis/) provides orbital-resolved band structure and DOS plots.

### orbvis (band-structure orbital visualization)

For periodic systems, [orbvis](https://pypi.org/project/orbvis/) provides orbital-resolved band structure and DOS plots.

```bash
pip install orbvis
```

### molecular-orbitals

[molecular-orbitals](https://pypi.org/project/molecular-orbitals/) provides a higher-level API for constructing and visualizing MOs from scratch using quantum chemistry methods.

```bash
pip install molecular-orbitals
```

### Hydrogen-like orbitals (analytical)

For pedagogical or analytical work, standard libraries are sufficient:

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm, genlaguerre, factorial

def hydrogen_wavefunction(n, l, m, r, theta, phi):
    """Hydrogen atom ψ_nlm in spherical coordinates (atomic units)."""
    a0 = 1.0  # Bohr radius
    rho = 2 * r / (n * a0)
    norm = np.sqrt(
        (2 / (n * a0))**3 * factorial(n - l - 1) /
        (2 * n * factorial(n + l)**3)
    )
    L = genlaguerre(n - l - 1, 2 * l + 1)(rho)
    Y = sph_harm(m, l, phi, theta)
    return norm * np.exp(-rho / 2) * rho**l * L * Y

# Plot ψ²_210 in the xz-plane
x = np.linspace(-20, 20, 400)
z = np.linspace(-20, 20, 400)
X, Z = np.meshgrid(x, z)
R = np.sqrt(X**2 + Z**2)
THETA = np.arctan2(np.sqrt(X**2), Z)

psi = hydrogen_wavefunction(2, 1, 0, R, THETA, np.zeros_like(R))
density = np.abs(psi)**2

plt.figure(figsize=(6, 6))
plt.contourf(X, Z, density, levels=80, cmap="inferno")
plt.colorbar(label="|ψ|²")
plt.title("2p_z orbital probability density")
plt.xlabel("x (a₀)")
plt.ylabel("z (a₀)")
plt.tight_layout()
plt.savefig("orbital_2pz.png", dpi=150)
plt.show()
```

---

## 2. Potential Energy Surfaces with Matplotlib and Plotly

### 2.1 Matplotlib – 2D contour and 3D surface

Matplotlib's `mpl_toolkits.mplot3d` handles static PES well and integrates directly with QC data.

**From raw energy grid (e.g. two dihedral angles):**
```python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

# Example: load energies from a 2D scan
# shape: (n_phi1, n_phi2)
phi1 = np.linspace(-180, 180, 36)
phi2 = np.linspace(-180, 180, 36)
PHI1, PHI2 = np.meshgrid(phi1, phi2)

# Replace with your actual energy array (in kcal/mol)
E = (np.cos(np.radians(PHI1)) + np.cos(np.radians(PHI2))
     + 0.5 * np.cos(np.radians(PHI1 - PHI2)))
E -= E.min()

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Contour map
cf = axes[0].contourf(PHI1, PHI2, E, levels=30, cmap="RdYlBu_r")
axes[0].contour(PHI1, PHI2, E, levels=30, colors="k", linewidths=0.4, alpha=0.5)
plt.colorbar(cf, ax=axes[0], label="ΔE (kcal/mol)")
axes[0].set_xlabel("φ₁ (°)")
axes[0].set_ylabel("φ₂ (°)")
axes[0].set_title("2D PES – Ramachandran-style")

# 3D surface
ax3d = fig.add_subplot(1, 2, 2, projection="3d")  # replaces axes[1]
ax3d.plot_surface(PHI1, PHI2, E, cmap="RdYlBu_r", edgecolor="none", alpha=0.85)
ax3d.set_xlabel("φ₁")
ax3d.set_ylabel("φ₂")
ax3d.set_zlabel("ΔE")
ax3d.set_title("3D PES")

plt.tight_layout()
plt.savefig("pes_matplotlib.png", dpi=150)
plt.show()
```

**Using autodE for automated QC-backed PES scans:**

[autodE](https://github.com/duartegroup/autode) wraps electronic structure codes and can generate 1D/2D surfaces directly.

```bash
pip install autode
```

```python
import autode as ade

ade.Config.n_cores = 8

# 2D relaxed scan over two breaking bonds (Diels-Alder example)
pes = ade.pes.RelaxedPESnD(
    species=ade.Molecule("reactant_complex.xyz"),
    rs={(0, 5): (3.0, 10), (3, 4): (3.0, 10)},  # atom pairs → (max_dist, n_steps)
)
pes.calculate(method=ade.methods.XTB())
pes.plot("DA_surface.png")          # matplotlib output
pes.save("DA_surface.npz")          # save for reuse

# Reload and plot with 4× spline interpolation
pes2 = ade.pes.RelaxedPESnD.from_file("DA_surface.npz")
pes2.plot("DA_surface_smooth.png", interp_factor=4)
```

### 2.2 Plotly – interactive 3D PES

Plotly excels at interactive, browser-renderable surfaces useful for exploration and presentations.

**Install:**
```bash
pip install plotly
```

**3D interactive surface (from numpy grid):**
```python
import numpy as np
import plotly.graph_objects as go

# Same grid as above; replace E with your data
phi1 = np.linspace(-180, 180, 72)
phi2 = np.linspace(-180, 180, 72)
PHI1, PHI2 = np.meshgrid(phi1, phi2)
E = (np.cos(np.radians(PHI1)) + np.cos(np.radians(PHI2))
     + 0.5 * np.cos(np.radians(PHI1 - PHI2)))
E -= E.min()

fig = go.Figure(data=[
    go.Surface(
        x=PHI1, y=PHI2, z=E,
        colorscale="RdBu",
        reversescale=True,
        colorbar=dict(title="ΔE (kcal/mol)"),
        contours=dict(
            z=dict(show=True, usecolormap=True, highlightcolor="white", project_z=True)
        ),
    )
])

fig.update_layout(
    title="Potential Energy Surface",
    scene=dict(
        xaxis_title="φ₁ (°)",
        yaxis_title="φ₂ (°)",
        zaxis_title="ΔE (kcal/mol)",
        camera=dict(eye=dict(x=1.5, y=1.5, z=1.0)),
    ),
    width=800, height=650,
)

fig.write_html("pes_interactive.html")   # open in browser
fig.show()
```

**Contour slice projection (useful for publication):**
```python
fig = go.Figure(data=[
    go.Surface(x=PHI1, y=PHI2, z=E, colorscale="Viridis",
               contours=dict(
                   x=dict(highlight=False),
                   y=dict(highlight=False),
                   z=dict(show=True, start=0, end=E.max(), size=0.5,
                          color="white", width=1),
               ))
])
fig.update_layout(scene_aspectmode="cube")
fig.show()
```

### 2.3 PESViewer – reaction network PES

For multi-well PES diagrams (wells, bimolecular products, TSs, barrierless reactions), use [PESViewer](https://github.com/rubenvdvijver/PESViewer):

```bash
pip install pesviewer
```

---

## 3. Reaction Paths

### 3.1 Kudi – IRC property extraction

[Kudi](https://github.com/stvogt/kudi) extracts chemically relevant properties (energy, atomic charges, bond lengths, etc.) along a Gaussian 16 IRC path and plots them cleanly.

**Install:**
```bash
pip install kudi
# or clone: git clone https://github.com/stvogt/kudi
```

**Basic IRC extraction and plotting:**
```python
from kudi import ReactionPath

# Point to a Gaussian IRC output file
rp = ReactionPath("reaction_irc.log")

# Access properties along the path
print(rp.energies)          # relative energies in kcal/mol
print(rp.irc_coordinate)    # IRC coordinate values

# Plot energy profile
rp.plot_energy(
    title="H-abstraction IRC",
    xlabel="IRC coordinate (amu^{1/2} bohr)",
    ylabel="ΔE (kcal/mol)",
    save="irc_energy.png",
)

# Plot any atomic property (e.g. NBO charge on atom 3) along the path
rp.plot_property(
    prop="nbo_charges",
    atom_index=3,
    title="NBO charge on O along IRC",
    save="irc_charge_O.png",
)
```

**Generate the IRC input if not yet computed:**
```bash
# Kudi ships a helper script for Gaussian 16
python make_irc_g16.py transition_state.fchk
```

### 3.2 Eyringpy – reaction force analysis

[Eyringpy](https://github.com/rperezbaamonde/eyringpy) extends IRC analysis with reaction force decomposition and activation strain.

```python
# After installing eyringpy
from eyringpy.irc import IRCAnalysis

irc = IRCAnalysis("irc_output.log", program="gaussian")
irc.reaction_force()           # F = -dE/dξ
irc.plot_force(save="rf.png")
irc.activation_strain()        # strain + interaction decomposition
```

### 3.3 Manual reaction path from energy + geometry arrays

When data comes from any source (ORCA, MOLPRO, etc.), build the path manually:

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# IRC coordinate and energies (populate from your parser)
xi = np.array([-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
E  = np.array([ 0.0,  5.2, 12.0, 22.5, 28.1, 19.3, 8.7, 2.1, 0.3])

cs = CubicSpline(xi, E)
xi_fine = np.linspace(xi[0], xi[-1], 500)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(xi_fine, cs(xi_fine), "b-", lw=2, label="Interpolated path")
ax.plot(xi, E, "ko", ms=6, label="QC points")
ax.axvline(0, color="gray", ls="--", lw=0.8, label="TS")
ax.fill_between(xi_fine, cs(xi_fine), alpha=0.1, color="blue")
ax.set_xlabel("IRC coordinate (amu$^{1/2}$ bohr)")
ax.set_ylabel("ΔE (kcal/mol)")
ax.set_title("Intrinsic Reaction Coordinate")
ax.legend()
plt.tight_layout()
plt.savefig("irc_manual.png", dpi=150)
plt.show()
```

### 3.4 GaussParse – parse Gaussian IRC output

[GaussParse](https://gaussparse.readthedocs.io/) provides utilities specifically for extracting IRC data from Gaussian files.

```bash
pip install gaussparse
```

---

## 4. Correlation Diagrams with PyBEST

[PyBEST](https://fizyka.umk.pl/~pybest/) (Pythonic Black-box Electronic Structure Tool) implements orbital entanglement analysis based on quantum information theory. It quantifies orbital correlations via the single-orbital entropy $s(1)_i$ and orbital-pair mutual information $I_{i|j}$, producing the correlation diagrams used in multireference analysis.

### Installation

PyBEST is distributed from the Toruń group's server. Check the [download page](https://fizyka.umk.pl/~pybest/pybest-v2.0.0/user_download_and_install.html) for the current distribution. Matplotlib ≥ 3.7 is required for visualization:

```bash
pip install "matplotlib~=3.7"
```

### Theory

The single-orbital entropy for orbital $i$ is:

$$s(1)_i = -\sum_{\alpha=1}^4 \omega_{\alpha,i} \ln \omega_{\alpha,i}$$

where $\omega_{\alpha,i}$ are eigenvalues of the one-orbital reduced density matrix. The orbital-pair mutual information is:

$$I_{i|j} = s(2)_{i,j} - s(1)_i - s(1)_j$$

Large $I_{i|j}$ indicates strong correlation between orbitals $i$ and $j$.

### Performing an orbital entanglement calculation (pCCD example)

```python
from pybest import context
from pybest.gbasis import (
    compute_eri, compute_kinetic, compute_nuclear,
    compute_nuclear_repulsion, compute_overlap, get_gobasis,
)
from pybest.geminals import ROOpCCD
from pybest.linalg import DenseLinalgFactory
from pybest.occ_model import AufbauOccModel
from pybest.orbital_entanglement import OrbitalEntanglementRpCCD
from pybest.wrappers import RHF

# --- Molecule and basis set ---
fn_xyz = context.get_fn("test/water.xyz")   # replace with your molecule
obasis  = get_gobasis("cc-pvdz", fn_xyz)

# --- Setup ---
lf        = DenseLinalgFactory(obasis.nbasis)
occ_model = AufbauOccModel(obasis, ncore=0)
orb_a     = lf.create_orbital(obasis.nbasis)
olp       = compute_overlap(obasis)

# --- Integrals ---
kin      = compute_kinetic(obasis)
ne       = compute_nuclear(obasis)
er       = compute_eri(obasis)
external = compute_nuclear_repulsion(obasis)

# --- HF reference ---
hf        = RHF(lf, occ_model)
hf_output = hf(kin, ne, er, external, olp, orb_a)

# --- Orbital-optimized pCCD ---
pccd        = ROOpCCD(lf, occ_model)
pccd_output = pccd(kin, ne, er, hf_output)

# --- Orbital entanglement ---
entanglement = OrbitalEntanglementRpCCD(lf, pccd_output)
entanglement()   # writes s1-pccd.dat and i12-pccd.dat
```

For pCCD-LCC corrections (higher accuracy):

```python
from pybest.cc import RpCCDLCCD, RpCCDLCCSD
from pybest.orbital_entanglement import OrbitalEntanglementRpCCDLCC

lccd        = RpCCDLCCD(lf, occ_model)
lccd_output = lccd(kin, ne, er, pccd_output, lambda_equations=True)

entanglement = OrbitalEntanglementRpCCDLCC(lf, lccd_output)
entanglement()   # writes s1-pccd-lcc.dat and i12-pccd-lcc.dat
```

### Generating the correlation diagram

PyBEST ships the `pybest-entanglement.py` script. After the calculation completes:

```bash
# Basic usage – uses default s1.dat and i12.dat
pybest-entanglement.py

# Specify method files and threshold
pybest-entanglement.py \
  --threshold 0.001 \
  --iname pybest-results/i12-pccd-lcc.dat \
  --sname pybest-results/s1-pccd-lcc.dat

# Restrict to at most 15 most-correlated orbitals
pybest-entanglement.py \
  --threshold 0.001 \
  --largest-i 15 \
  --iname pybest-results/i12-pccd-lcc.dat \
  --sname pybest-results/s1-pccd-lcc.dat
```

**Output:** `pybest-results/s1.pdf` (single-orbital entropy bar chart) and `pybest-results/i12.pdf` (mutual information network diagram). If MO pictures are saved as `mo_[index].png`, the script automatically arranges them around the network.

**Key script arguments:**

| Argument | Type | Description |
|---|---|---|
| `--threshold` | float | Minimum $I_{i|j}$ to display (e.g. `0.001`) |
| `--iname` | str | Path to mutual information `.dat` file |
| `--sname` | str | Path to single-orbital entropy `.dat` file |
| `--indices` | int… | Restrict to specific orbital indices |
| `--order` | int… | Custom orbital ordering in the plot |
| `--zoom` | float | Scale factor for MO pictures |
| `--largest-i` | int | Keep only N most correlated orbitals |

---

## Quick Reference

| Task | Tool | Output format |
|---|---|---|
| Molecular orbital isosurfaces (interactive/static/batch) | orbital-viz | HTML / png (Plotly) |
| Molecular orbital isosurfaces (static/batch) | ORBKIT + PyVista | PNG / HDF5 |
| Analytical hydrogen orbitals | NumPy + Matplotlib | PNG |
| 2D/3D PES from grid data | Matplotlib `mplot3d` | PNG |
| Automated QC PES scans | autodE | PNG / NPZ |
| Interactive PES | Plotly `go.Surface` | HTML |
| Multi-well reaction network PES | PESViewer | SVG / PNG |
| IRC property extraction | Kudi | PNG |
| Reaction force analysis | Eyringpy | PNG |
| Orbital correlation / entanglement | PyBEST | PDF |

---


## Interpretation hints for orbitals

- **Bonding MO**: large amplitude _between_ atoms, lobes in phase → constructive interference
- **Antibonding MO**: nodal plane _between_ bonded atoms, lobes out of phase
- **Non-bonding MO**: amplitude localised on one atom / lone pair, minimal overlap with neighbours
- **π vs σ**: π MOs have a nodal plane containing the molecular axis; σ MOs are cylindrically symmetric along the bond
- Isovalue 0.05 captures ~85–90 % of the electron density in most MOs; 0.02 reveals diffuse regions


## References

- orbital-viz: [pypi.org/project/orbital-viz](https://pypi.org/project/orbital-viz/) | [github.com/AleksandraLeszczyk/OrbitalViz](https://github.com/AleksandraLeszczyk/OrbitalViz)
- ORBKIT: [orbkit.github.io](https://orbkit.github.io/)
- autodE: [duartegroup.github.io/autodE](https://duartegroup.github.io/autodE)
- PESViewer: [github.com/rubenvdvijver/PESViewer](https://github.com/rubenvdvijver/PESViewer)
- Kudi: [github.com/stvogt/kudi](https://github.com/stvogt/kudi) | [doi:10.1007/s00894-016-2983-3](https://doi.org/10.1007/s00894-016-2983-3)
- Eyringpy: [github.com/rperezbaamonde/eyringpy](https://github.com/rperezbaamonde/eyringpy)
- PyBEST: [fizyka.umk.pl/~pybest](https://fizyka.umk.pl/~pybest/pybest-v2.0.0/index.html) | [doi:10.1016/j.cpc.2021.107933](https://doi.org/10.1016/j.cpc.2021.107933)
- Boguslawski & Tecmer, orbital entanglement: [doi:10.1002/qua.24832](https://doi.org/10.1002/qua.24832)
