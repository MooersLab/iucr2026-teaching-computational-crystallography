# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Demo 1. DIALS in Jupyter
#
# **Companion to slide 7 of *Teaching Crystallographic Computing with Computational Notebooks*.**
# 27th Congress and General Assembly of the IUCr, Calgary, Alberta, 2026 August 17.
#
# Blaine Mooers, PhD. Department of Biochemistry and Physiology, University of Oklahoma Health Campus.
# blaine-mooers@ou.edu
#
# ## What this demonstration covers
#
# 1. Importing diffraction images.
# 2. Spot finding and indexing.
# 3. Refinement and integration.
# 4. Visualizing reciprocal space.
#
# ## The idea behind the demonstration
#
# DIALS is a framework rather than a monolithic integration program, so every step of
# data reduction is a separate dispatcher that reads and writes files you can inspect.
# That design is a gift to teaching. Students watch one operation at a time, look at what
# it produced, and only then move on. The pipeline stops being a single opaque button.
#
# Each dispatcher also has a Python counterpart. The notebook uses the shell form for the
# workflow and the Python form for the inspection, so students see both faces of the same
# library.
#
# ## Data
#
# Zenodo record [3752540](https://doi.org/10.5281/zenodo.3752540) holds 686 CBF frames of
# SARS-CoV-2 main protease bound to 2-methyl-1-tetralone. The frames were collected on
# beamline P11 at PETRA III with 0.2 degrees per frame, a 40 ms exposure, a wavelength of
# 1.033214 angstroms, and the detector at 200 mm. The deposited structure is PDB entry 6YNQ.
#
# Each frame is 6.2 MB, so the full sweep is 6.2 GB. We take a wedge instead.

# %% [markdown]
# ## How to run this in front of a class
#
# The conda install of DIALS takes ten to fifteen minutes on a cold Colab virtual machine.
# Run this notebook once the night before, save a copy to your Google Drive, and open that
# copy in class. Then only the workflow cells need to run live.
#
# Cells marked **PREBAKE** are the slow ones.
#
# If a cell stalls, use *Runtime > Restart runtime* and rerun from the top of the section.

# %% [markdown]
# ## 1. Where are we?
#
# Colab hands you a fresh Linux virtual machine for every session. The machine disappears
# when you disconnect, so anything worth keeping belongs on Google Drive or on your laptop.

# %%
# Silence the pydevd file-validation warning that Colab prints on import.
# This must be set before anything else imports.
import os

os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"

import platform
import shutil
import sys

print("Python   :", sys.version.split()[0])
print("Platform :", platform.platform())
try:
    total, used, free = shutil.disk_usage("/content")
    root = "/content"
except FileNotFoundError:
    total, used, free = shutil.disk_usage(".")
    root = os.getcwd()
print(f"Disk on {root}: {free / 1e9:.1f} GB free of {total / 1e9:.1f} GB")

IN_COLAB = "google.colab" in sys.modules or os.path.isdir("/content")
print("Running in Colab:", IN_COLAB)

# %%
# Make a clean working directory and move into it.
import pathlib

WORK = pathlib.Path("/content/mpro_dials" if IN_COLAB else "./mpro_dials").resolve()
WORK.mkdir(parents=True, exist_ok=True)
os.chdir(WORK)
print("Working in", WORK)

# %% [markdown]
# ## 2. Install DIALS
#
# **PREBAKE.** DIALS ships on conda-forge together with cctbx, dxtbx, and scitbx. The
# cleanest route onto Colab is [condacolab](https://github.com/conda-incubator/condacolab),
# which swaps the system Python for a Mambaforge install and then restarts the kernel.
#
# The kernel restart in the next cell is expected. Do not panic and do not rerun the cell.

# %%
# PREBAKE. Installs condacolab and restarts the kernel.
if IN_COLAB:
    # !pip -q install condacolab
    import condacolab

    condacolab.install()  # The kernel restarts here. This is normal.
else:
    print("Not in Colab. Install DIALS with: conda create -n dials -c conda-forge dials")

# %%
# PREBAKE. Run this only after the kernel has restarted.
#
# condacolab writes a pinned file that locks python and python_abi.
# Conda-forge resolves the DIALS stack against a different minor version,
# so wipe the pin and let the solver move freely. This is safe inside a
# throwaway Colab virtual machine.
# !rm -f /usr/local/conda-meta/pinned
# !mamba install -y -q -c conda-forge dials gemmi reciprocalspaceship matplotlib pandas

# %%
# Smoke test. dials.version prints the build identifier.
# !dials.version
# !which dials.import dials.find_spots dials.index dials.refine dials.integrate dials.scale

# %% [markdown]
# ### Re-establish the session after the restart
#
# The kernel restart wiped every variable. This cell rebuilds the three names the rest of
# the notebook depends on. Run it once after the DIALS install finishes.

# %%
import os
import pathlib
import sys

os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"
IN_COLAB = "google.colab" in sys.modules or os.path.isdir("/content")
WORK = pathlib.Path("/content/mpro_dials" if IN_COLAB else "./mpro_dials").resolve()
WORK.mkdir(parents=True, exist_ok=True)
os.chdir(WORK)
print("Working in", WORK)

# %% [markdown]
# ## 3. Fetch a wedge of frames from Zenodo
#
# Sixty frames span twelve degrees of rotation. That is enough to find spots, index, and
# integrate a partial dataset, and it downloads in a few minutes rather than an hour.
# Raise `N_FRAMES` when you want better completeness and you have the time to wait.

# %%
import concurrent.futures
import urllib.request

RECORD = "3752540"
STEM = "l6p17_09_001"
N_FRAMES = 60  # 60 frames x 0.2 deg = 12 degrees of rotation.
FRAME_DIR = WORK / "frames"
FRAME_DIR.mkdir(parents=True, exist_ok=True)


def frame_url(n: int) -> str:
    """Return the direct download URL for frame n of the Zenodo record."""
    return f"https://zenodo.org/records/{RECORD}/files/{STEM}_{n:05d}.cbf?download=1"


import time


def fetch_frame(n: int, attempts: int = 4) -> str:
    """Download frame n unless it is already on disk. Retry on transient errors.

    urllib.request.urlretrieve writes chunks straight to disk, so a mid-stream
    drop leaves a truncated file. We delete the partial file before each retry
    so the next attempt starts clean.
    """
    dest = FRAME_DIR / f"{STEM}_{n:05d}.cbf"
    # A cached frame is accepted only if it is within 0.5 percent of the
    # largest frame already on disk. This catches truncated downloads that
    # sneak past a fixed byte threshold, because Zenodo can drop a
    # connection after megabytes of data have already streamed to disk.
    peers = [p.stat().st_size for p in FRAME_DIR.glob("*.cbf")]
    peak = max(peers) if peers else 0
    if dest.exists() and dest.stat().st_size >= 0.995 * peak and peak > 0:
        return f"cached  {dest.name}"
    if dest.exists():
        dest.unlink()  # remove the partial file before retrying
    last_err = None
    for attempt in range(1, attempts + 1):
        try:
            urllib.request.urlretrieve(frame_url(n), dest)
            return f"fetched {dest.name} ({dest.stat().st_size / 1e6:.1f} MB)"
        except Exception as exc:  # ContentTooShortError, URLError, timeout
            last_err = exc
            if dest.exists():
                dest.unlink()
            time.sleep(2 ** attempt)  # 2, 4, 8, 16 s backoff
    return f"FAILED  {dest.name} after {attempts} attempts: {last_err}"


# Four workers is polite to Zenodo. Eight is where the throttling starts to bite.
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    for msg in pool.map(fetch_frame, range(1, N_FRAMES + 1)):
        print(msg)

# %%
# Confirm the wedge is complete and readable.
frames = sorted(FRAME_DIR.glob(f"{STEM}_*.cbf"))
total_mb = sum(f.stat().st_size for f in frames) / 1e6
print(f"{len(frames)} frames, {total_mb:.0f} MB")
print("First :", frames[0].name)
print("Last  :", frames[-1].name)

# %% [markdown]
# ### Look at a frame before processing it
#
# Students should see the diffraction pattern before any program touches it. A CBF frame is
# a text header followed by a compressed image, and `dxtbx` reads both.

# %%
import matplotlib.pyplot as plt
import numpy as np
from dxtbx.model.experiment_list import ExperimentListFactory

expts_one = ExperimentListFactory.from_filenames([str(frames[0])])
imageset = expts_one.imagesets()[0]
raw = imageset.get_raw_data(0)[0].as_numpy_array()

fig, ax = plt.subplots(figsize=(7, 7))
ax.imshow(np.clip(raw, 0, 60), cmap="Greys", origin="lower")
ax.set_title(f"{frames[0].name}\nintensities clipped at 60 counts")
ax.set_xlabel("fast axis, pixels")
ax.set_ylabel("slow axis, pixels")
plt.tight_layout()
plt.show()

# %%
# The same frame described by its models rather than its pixels.
beam = imageset.get_beam()
det = imageset.get_detector()[0]
gon = imageset.get_goniometer()
scan = imageset.get_scan()

print("Wavelength      :", round(beam.get_wavelength(), 6), "angstroms")
print("Detector origin :", tuple(round(v, 3) for v in det.get_origin()), "mm")
print("Pixel size      :", det.get_pixel_size(), "mm")
print("Image size      :", det.get_image_size(), "pixels")
print("Distance        :", round(det.get_distance(), 2), "mm")
print("Rotation axis   :", gon.get_rotation_axis())
print("Oscillation     :", scan.get_oscillation(), "degrees")

# %% [markdown]
# ## 4. `dials.import`
#
# The first dispatcher reads the headers of every frame and writes `imported.expt`, a JSON
# file holding the beam, detector, goniometer, and scan models. No pixel data is copied.

# %%
# !dials.import frames/{STEM}_*.cbf output.experiments=imported.expt

# %%
# imported.expt is plain JSON. Show students that nothing is hidden.
import json

with open("imported.expt") as fh:
    doc = json.load(fh)
print("Top-level keys :", list(doc.keys()))
print()
print(json.dumps(doc["beam"], indent=2)[:600])

# %% [markdown]
# ## 5. `dials.find_spots`
#
# Spot finding applies a dispersion threshold to every frame and records the position and
# intensity of each strong reflection. The result is `strong.refl`, a reflection table.
#
# This is the first place where a parameter change produces a visible difference, so it is
# the first good place to stop and let a student change something.

# %%
# !dials.find_spots imported.expt output.reflections=strong.refl nproc=2

# %%
# A reflection table behaves like a column store. Print its shape and columns.
from dials.array_family import flex

strong = flex.reflection_table.from_file("strong.refl")
print("Strong spots :", strong.size())
print("Columns      :")
for name in sorted(strong.keys()):
    print("   ", name)

# %%
# Where did the spots land on the detector, and how do they spread through the sweep?
xyz = strong["xyzobs.px.value"].as_double().as_numpy_array().reshape(-1, 3)
intensity = strong["intensity.sum.value"].as_numpy_array()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

sc = ax1.scatter(xyz[:, 0], xyz[:, 1], c=np.log10(np.clip(intensity, 1, None)), s=3, cmap="viridis")
ax1.set_aspect("equal")
ax1.set_xlabel("fast axis, pixels")
ax1.set_ylabel("slow axis, pixels")
ax1.set_title(f"{strong.size()} strong spots on the detector")
fig.colorbar(sc, ax=ax1, label="log10 summed intensity")

ax2.hist(xyz[:, 2], bins=N_FRAMES)
ax2.set_xlabel("image number")
ax2.set_ylabel("strong spots")
ax2.set_title("Spots per frame")
plt.tight_layout()
plt.show()

# %% [markdown]
# **Try this in class.** Rerun spot finding with `spotfinder.threshold.dispersion.sigma_strong=5`
# and compare the spot count. A higher threshold finds fewer, cleaner spots. Ask the students
# which direction helps indexing and which direction helps integration.

# %% [markdown]
# ## 6. Visualizing reciprocal space
#
# Spot finding gives positions on a detector. Mapping those positions into reciprocal space
# turns them into a lattice, and the lattice is the thing indexing actually works on. Showing
# the lattice before indexing runs makes the next step feel inevitable rather than magical.

# %%
from dxtbx.model.experiment_list import ExperimentListFactory

expts = ExperimentListFactory.from_json_file("imported.expt", check_format=False)

mapped = flex.reflection_table.from_file("strong.refl")
mapped.centroid_px_to_mm(expts)
mapped.map_centroids_to_reciprocal_space(expts)
rlp = mapped["rlp"].as_double().as_numpy_array().reshape(-1, 3)
print("Reciprocal lattice points :", rlp.shape)

# %%
# Three orthogonal slabs through the reciprocal lattice.
fig, axes = plt.subplots(1, 3, figsize=(16, 5.4))
slab = 0.004  # thickness of the slab in reciprocal angstroms
pairs = [(0, 1, 2, "$x^*$", "$y^*$"), (0, 2, 1, "$x^*$", "$z^*$"), (1, 2, 0, "$y^*$", "$z^*$")]

for ax, (i, j, k, xl, yl) in zip(axes, pairs):
    keep = np.abs(rlp[:, k]) < slab
    ax.scatter(rlp[keep, i], rlp[keep, j], s=2, c="k")
    ax.set_aspect("equal")
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(f"slab through {yl.strip('$')}, {keep.sum()} points")

fig.suptitle("Reciprocal lattice from the strong spots, before indexing")
plt.tight_layout()
plt.show()

# %% [markdown]
# The rows of dots are the reciprocal lattice. Indexing is the search for the three vectors
# that generate them.

# %% [markdown]
# ## 7. `dials.index`
#
# Indexing finds the orientation matrix and the unit cell. The default one-dimensional
# Fourier transform method searches for periodic spacing along many candidate directions.

# %%
# !dials.index imported.expt strong.refl output.experiments=indexed.expt output.reflections=indexed.refl

# %%
# Pull the crystal model out of the experiment list rather than reading it off the log.
indexed_expts = ExperimentListFactory.from_json_file("indexed.expt", check_format=False)
xtal = indexed_expts[0].crystal

print("Space group   :", xtal.get_space_group().info())
uc = xtal.get_unit_cell().parameters()
print("Unit cell     : a={:.2f} b={:.2f} c={:.2f} alpha={:.1f} beta={:.1f} gamma={:.1f}".format(*uc))
print("Volume        : {:.0f} cubic angstroms".format(xtal.get_unit_cell().volume()))
print()
print("A matrix      :")
print(np.array(xtal.get_A()).reshape(3, 3))

# %%
# How many of the strong spots picked up an index?
indexed = flex.reflection_table.from_file("indexed.refl")
flags = indexed.get_flags(indexed.flags.indexed)
print(f"{flags.count(True)} of {indexed.size()} reflections were indexed "
      f"({100 * flags.count(True) / indexed.size():.1f} percent)")

# %% [markdown]
# **Discussion point.** The deposited structure 6YNQ is in space group C2. `dials.index`
# reports the lattice it can defend from the geometry alone, which is usually the primitive
# triclinic setting or the highest metric symmetry consistent with the cell. Symmetry is
# determined later, from the intensities, by `dials.symmetry`. Separating the two steps is
# worth a minute of class time, because students routinely conflate them.

# %% [markdown]
# ## 8. `dials.refine`
#
# Refinement polishes the experimental geometry against every indexed reflection at once.
# The detector position, the beam direction, and the crystal orientation all move.

# %%
# !dials.refine indexed.expt indexed.refl output.experiments=refined.expt output.reflections=refined.refl

# %%
# Compare the cell before and after refinement.
refined_expts = ExperimentListFactory.from_json_file("refined.expt", check_format=False)
before = indexed_expts[0].crystal.get_unit_cell().parameters()
after = refined_expts[0].crystal.get_unit_cell().parameters()

print(f"{'':10s}{'indexed':>12s}{'refined':>12s}{'change':>12s}")
for label, b, a in zip(["a", "b", "c", "alpha", "beta", "gamma"], before, after):
    print(f"{label:10s}{b:12.4f}{a:12.4f}{a - b:12.4f}")

# %% [markdown]
# ## 9. `dials.integrate` and `dials.scale`
#
# Integration measures the intensity of every predicted reflection, including the weak ones
# that spot finding skipped. Scaling puts the measurements on a common scale and merges them.
#
# Integration is the slowest workflow step. On two Colab cores a twelve-degree wedge takes a
# few minutes.

# %%
# !dials.integrate refined.expt refined.refl output.experiments=integrated.expt output.reflections=integrated.refl nproc=2

# %%
# !dials.symmetry integrated.expt integrated.refl output.experiments=symmetrized.expt output.reflections=symmetrized.refl

# %%
# !dials.scale symmetrized.expt symmetrized.refl output.html=scaling.html output.unmerged_mtz=scaled_unmerged.mtz output.merged_mtz=scaled.mtz

# %% [markdown]
# ## 10. Read the merged data back without leaving Python
#
# The point of ending here is that the MTZ file is not a terminus. It is an input to the
# next notebook, and it is a pandas data frame away from any analysis a student wants to try.

# %%
import gemmi

mtz = gemmi.read_mtz_file("scaled.mtz")
print("Space group :", mtz.spacegroup.hm)
print("Cell        :", mtz.cell)
print("Resolution  : {:.2f} to {:.2f} angstroms".format(mtz.resolution_low(), mtz.resolution_high()))
print("Columns     :")
for col in mtz.columns:
    print(f"    {col.label:12s} {col.type}")

# %%
# Signal-to-noise against resolution, computed in twenty shells.
import reciprocalspaceship as rs

ds = rs.read_mtz("scaled.mtz").compute_dHKL()
ds["invd2"] = 1.0 / ds["dHKL"] ** 2
ds["shell"] = np.digitize(ds["invd2"], np.linspace(ds["invd2"].min(), ds["invd2"].max(), 21))

icol = "IMEAN" if "IMEAN" in ds.columns else "I"
scol = "SIGIMEAN" if "SIGIMEAN" in ds.columns else "SIGI"
agg = ds.groupby("shell").agg(I=(icol, "mean"), SIG=(scol, "mean"), d=("dHKL", "mean"))
agg["IsigI"] = agg["I"] / agg["SIG"]

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(1.0 / agg["d"] ** 2, agg["IsigI"], marker="o")
ax.axhline(2.0, ls="--", c="grey", label="I/sigma(I) = 2")
ax.set_xlabel(r"$1/d^2$, inverse square angstroms")
ax.set_ylabel(r"$\langle I \rangle / \langle \sigma(I) \rangle$")
ax.set_title("Data quality against resolution")
ax.legend()

# Label the top axis in angstroms, because that is how crystallographers think.
sec = ax.secondary_xaxis("top", functions=(lambda s: np.sqrt(1 / np.maximum(s, 1e-9)),
                                           lambda d: 1 / np.maximum(d, 1e-9) ** 2))
sec.set_xlabel("resolution, angstroms")
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 11. The DIALS report
#
# `dials.report` collects every diagnostic plot into one HTML file. Download it and open it
# in a browser tab beside the notebook.

# %%
# !dials.report symmetrized.expt symmetrized.refl output.html=dials.report.html

# %%
if IN_COLAB:
    from google.colab import files

    files.download("dials.report.html")
    files.download("scaled.mtz")
else:
    print("Report written to", (WORK / "dials.report.html").resolve())

# %% [markdown]
# ## Exercises
#
# 1. Rerun spot finding with `sigma_strong=5` and again with `sigma_strong=2`. Plot the spot
#    count against the threshold. Which setting gives the best indexing rate?
# 2. Cut the wedge in half by importing only the first thirty frames. Does indexing still
#    succeed? What happens to completeness after scaling?
# 3. Force a wrong lattice with `dials.index space_group=P1` and compare the refined cell
#    against the default run.
# 4. Use `mapped["rlp"]` to compute the resolution of every strong spot, then plot the
#    distribution. Compare it against the `resolution: 1.304257A` claim in the Zenodo header.
#
# ## What to carry forward
#
# `scaled.mtz` from this notebook is the input to Demo 3, the molecular replacement workflow.
# Save it to Google Drive if you plan to run these demonstrations in sequence.

# %%
# Keep the merged data after the virtual machine shuts down.
# Set SAVE_TO_DRIVE to True the first time you run this in class.
SAVE_TO_DRIVE = False

if IN_COLAB and SAVE_TO_DRIVE:
    import shutil

    from google.colab import drive

    drive.mount("/content/drive")
    dest = pathlib.Path("/content/drive/MyDrive/iucr26_demos")
    dest.mkdir(parents=True, exist_ok=True)
    for name in ["scaled.mtz", "scaled_unmerged.mtz"]:
        shutil.copy(name, dest / name)
        print("Copied", name, "to", dest)

# %% [markdown]
# ## Sources
#
# * DIALS documentation. https://dials.github.io
# * Zenodo record 3752540. https://doi.org/10.5281/zenodo.3752540
# * PDB entry 6YNQ. https://www.rcsb.org/structure/6YNQ
