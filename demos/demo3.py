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
# # Demo 3. Molecular replacement workflow
#
# **Companion to slide 10 of *Teaching Crystallographic Computing with Computational Notebooks*.**
# 27th Congress and General Assembly of the IUCr, Calgary, Alberta, 2026 August 17.
#
# Blaine Mooers, PhD. Department of Biochemistry and Physiology, University of Oklahoma Health Campus.
# blaine-mooers@ou.edu
#
# ## What this demonstration covers
#
# 1. Preparing search models with CHAINSAW.
# 2. Running MOLREP or PHASER.
# 3. Evaluating MR solutions.
# 4. Initial refinement cycles.
#
# ## An honest note about what runs where
#
# The full CCP4 suite is six to eight gigabytes and its graphical front end does not run
# inside a notebook. Phenix requires a personalized installer tarball. Neither fits a free
# Colab session, so this notebook splits into two kinds of cell.
#
# * **Shown, not run.** The MOLREP, PHASER, and CHAINSAW command lines appear exactly as you
#   would type them on a workstation. Copy them out and run them there.
# * **Run here.** Model preparation, solution evaluation, and rigid-body refinement all run
#   in Colab on nothing heavier than `gemmi`, `numpy`, and `scipy`.
#
# The part that students find hardest is not launching Phaser. Launching Phaser is one
# command. The hard part is answering *is this solution right?*, and that part runs here in
# full.
#
# ## Data
#
# * **Target.** SARS-CoV-2 main protease bound to 2-methyl-1-tetralone, PDB entry
#   [6YNQ](https://www.rcsb.org/structure/6YNQ). Demo 1 processes the raw frames for this
#   structure. If you have not run Demo 1, this notebook downloads the deposited structure
#   factors instead.
# * **Search model.** Apo main protease, PDB entry
#   [6Y2E](https://www.rcsb.org/structure/6Y2E).

# %%
import os

os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"

import pathlib
import sys

IN_COLAB = "google.colab" in sys.modules or os.path.isdir("/content")
WORK = pathlib.Path("/content/mpro_mr" if IN_COLAB else "./mpro_mr").resolve()
WORK.mkdir(parents=True, exist_ok=True)
os.chdir(WORK)
print("In Colab   :", IN_COLAB)
print("Working in :", WORK.resolve())

# %%
# Everything this notebook needs installs in seconds.
# !pip -q install gemmi reciprocalspaceship scipy matplotlib pandas

# %%
import gemmi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import reciprocalspaceship as rs
import scipy.optimize

print("gemmi :", gemmi.__version__)
print("rs    :", rs.__version__)

# %% [markdown]
# ## 1. Get the observed amplitudes
#
# Two routes. Prefer the MTZ that Demo 1 produced, because processing your own data and then
# phasing it is the point of the exercise. Fall back to the deposited structure factors when
# Demo 1 has not been run.

# %%
import subprocess
import urllib.request

# DATA_SOURCE picks which reflections feed molecular replacement.
#   "auto"  : report what is available, then default to "pdb" because a
#             12-degree wedge from Demo 1 is far too incomplete for MR.
#   "demo1" : force the merged data from Demo 1. Will almost certainly fail
#             the rotation function unless you re-ran Demo 1 with a much
#             larger sweep and got completeness above about 80 percent.
#   "pdb"   : force the deposited structure factors for 6YNQ.
DATA_SOURCE = "auto"

demo1_mtz = pathlib.Path(
    "/content/mpro_dials/scaled.mtz" if IN_COLAB else "../mpro_dials/scaled.mtz"
)

# Report on the Demo 1 file if it exists, so the user knows what they would
# be picking if they overrode DATA_SOURCE.
n_demo1_refl = 0
if demo1_mtz.exists():
    demo1_probe = gemmi.read_mtz_file(str(demo1_mtz))
    n_demo1_refl = demo1_probe.nreflections
    d_min = demo1_probe.resolution_high()
    print(f"Demo 1 MTZ found : {n_demo1_refl} reflections at d_min = {d_min:.2f} A")
    print("A full Mpro dataset at that resolution has roughly 65,000 unique reflections.")
    print("A 12-degree wedge yields about 5 to 8 percent completeness, well below the")
    print("60 to 80 percent that molecular replacement needs.")
else:
    print("Demo 1 MTZ not found on disk.")

# Auto mode picks "demo1" only if the merged file plausibly reaches 60 percent
# completeness for a monoclinic cell at 1.3 A. Otherwise it picks "pdb".
if DATA_SOURCE == "auto":
    choice = "demo1" if n_demo1_refl >= 40_000 else "pdb"
    print(f"\nAuto-selected data source: {choice.upper()}")
    print("(set DATA_SOURCE at the top of this cell to override)")
else:
    choice = DATA_SOURCE
    print(f"\nDATA_SOURCE = {choice.upper()}")

OBS_MTZ = WORK / "obs.mtz"
if choice == "demo1":
    OBS_MTZ.write_bytes(demo1_mtz.read_bytes())
    print("Using the merged data from Demo 1.")
elif choice == "pdb":
    sf_cif = WORK / "6YNQ-sf.cif"
    if not sf_cif.exists():
        urllib.request.urlretrieve("https://files.rcsb.org/download/6YNQ-sf.cif", sf_cif)
    # gemmi converts the deposited structure-factor CIF into an MTZ.
    subprocess.run(["gemmi", "cif2mtz", str(sf_cif), str(OBS_MTZ)], check=True)
    print("Using the deposited structure factors for 6YNQ.")
else:
    raise ValueError(f"DATA_SOURCE must be 'auto', 'demo1', or 'pdb'; got {DATA_SOURCE!r}")

mtz = gemmi.read_mtz_file(str(OBS_MTZ))
print("Space group :", mtz.spacegroup.hm)
print("Cell        :", mtz.cell)
print("Resolution  : {:.2f} to {:.2f} angstroms".format(mtz.resolution_low(), mtz.resolution_high()))
print("Columns     :", [f"{c.label}({c.type})" for c in mtz.columns])

# %%
# Reduce whatever we were given to a single amplitude column named FOBS.
ds = rs.read_mtz(str(OBS_MTZ))

if "FP" in ds.columns:
    ds["FOBS"], ds["SIGFOBS"] = ds["FP"], ds["SIGFP"]
elif "F-obs-filtered" in ds.columns:
    ds["FOBS"], ds["SIGFOBS"] = ds["F-obs-filtered"], ds["SIGF-obs-filtered"]
elif "F" in ds.columns:
    ds["FOBS"], ds["SIGFOBS"] = ds["F"], ds["SIGF"]
else:
    # Only intensities were merged. Take the square root as a teaching-grade
    # substitute for a proper French and Wilson treatment.
    icol = "IMEAN" if "IMEAN" in ds.columns else "I"
    scol = "SIGIMEAN" if "SIGIMEAN" in ds.columns else "SIGI"
    keep = ds[icol] > 0
    ds = ds.loc[keep].copy()
    ds["FOBS"] = np.sqrt(ds[icol].to_numpy())
    ds["SIGFOBS"] = ds[scol].to_numpy() / (2 * ds["FOBS"].to_numpy())
    print("Converted intensities to amplitudes by taking square roots.")
    print("On real data run phenix.french_wilson or ctruncate instead.")

ds = ds[["FOBS", "SIGFOBS"]].dropna()
print(f"{len(ds)} reflections carried forward")
ds.head()

# %% [markdown]
# ## 2. Prepare the search model
#
# ### What CHAINSAW does, and why
#
# A search model that is too complete is worse than one that is trimmed. Side chains that
# differ between the model and the target contribute Fcalc in the wrong places, which lowers
# the signal in the rotation and translation functions. CHAINSAW takes an alignment between
# the search model sequence and the target sequence and does one of three things to each
# residue.
#
# * **Identical residue.** Keep every atom.
# * **Different residue.** Truncate to the last common atom, usually CB or CG.
# * **No equivalent.** Delete the residue.
#
# On a workstation with CCP4 installed, the run is this.
#
# ```bash
# chainsaw XYZIN 6Y2E.pdb XYZOUT 6Y2E_chainsaw.pdb ALIGNIN mpro.aln <<EOF
# MODE MIXS
# END
# EOF
# ```
#
# `MODE MIXS` keeps identical side chains and truncates the rest. `MODE MIXA` truncates every
# non-identical residue to alanine. `MODE MIXC` truncates everything to a poly-serine.
#
# Here the search model and the target are the same protein, so the alignment step is
# unnecessary and the useful trimming is different. We remove waters, ligands, hydrogens, and
# alternate conformers, then set every occupancy to one. Doing that with `gemmi` in eight
# lines shows students exactly what a model-preparation program is doing.

# %%
for pdb_id in ["6Y2E", "6YNQ"]:
    dest = WORK / f"{pdb_id}.cif"
    if not dest.exists():
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.cif", dest)
    print(f"{pdb_id}.cif  {dest.stat().st_size / 1e3:.0f} kB")

# %%
def prepare_search_model(path_in, path_out):
    """Strip a deposited structure down to a molecular replacement search model."""
    st = gemmi.read_structure(str(path_in))
    st.setup_entities()
    before = st[0].count_atom_sites()

    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    st.remove_empty_chains()

    for model in st:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.occ = 1.0

    after = st[0].count_atom_sites()
    st.setup_entities()
    st.write_pdb(str(path_out))
    print(f"{path_in.name}: {before} atoms in, {after} atoms out, {before - after} removed")
    return st


search = prepare_search_model(WORK / "6Y2E.cif", WORK / "6Y2E_search.pdb")
target = gemmi.read_structure(str(WORK / "6YNQ.cif"))
target.setup_entities()

# %%
# A poly-alanine variant, the CHAINSAW MIXA equivalent. Useful when the search model
# is a distant homologue rather than the same protein.
#
# gemmi ships `Residue.trim_to_alanine()`, so the whole operation is four lines. The
# commented-out block below does the same thing by hand, which is worth showing once
# so students see that "truncate the side chain" means "delete atoms from a list".


def truncate_to_polyalanine(st_in, path_out):
    st = st_in.clone()
    for model in st:
        for chain in model:
            for residue in chain:
                residue.trim_to_alanine()
                # Equivalent by hand:
                # keep = [a for a in residue if a.name in {"N", "CA", "C", "O", "CB"}]
                # while len(residue):
                #     del residue[0]
                # for atom in keep:
                #     residue.add_atom(atom)
    st.write_pdb(str(path_out))
    print("Poly-alanine model :", st[0].count_atom_sites(), "atoms")
    return st


poly_ala = truncate_to_polyalanine(search, WORK / "6Y2E_polyala.pdb")

# %% [markdown]
# ## 3. The molecular replacement search itself
#
# Copy one of these onto a machine that has the software. Nothing in this section runs in
# Colab, and pretending otherwise would be dishonest.
#
# ### PHASER, through Phenix
#
# ```bash
# phenix.phaser \
#     hklin=obs.mtz \
#     labin="F=FOBS SIGF=SIGFOBS" \
#     model=6Y2E_search.pdb \
#     model.identity=99 \
#     composition.chain.sequence=mpro.seq \
#     composition.chain.copies=1 \
#     sgalternative=all \
#     output.directory=phaser_out
# ```
#
# ### PHASER, through CCP4
#
# ```bash
# phaser <<EOF
# MODE MR_AUTO
# HKLIN obs.mtz
# LABIN F=FOBS SIGF=SIGFOBS
# ENSEMBLE mpro PDBFILE 6Y2E_search.pdb IDENTITY 0.99
# COMPOSITION PROTEIN SEQUENCE mpro.seq NUMBER 1
# SEARCH ENSEMBLE mpro NUMBER 1
# ROOT phaser_out
# EOF
# ```
#
# ### MOLREP
#
# ```bash
# molrep HKLIN obs.mtz MODEL 6Y2E_search.pdb <<EOF
# LABIN F=FOBS SIGF=SIGFOBS
# NMON 1
# END
# EOF
# ```
#
# ### Without installing anything
#
# [CCP4 Cloud](https://cloud.ccp4.ac.uk) runs the same programs in a browser tab. Move files
# between Colab and CCP4 Cloud through Google Drive. For a class of twenty this is often the
# path of least resistance, because nobody has to install anything.
#
# ### Driving any of them from a notebook
#
# The `run` helper from Demo 2 is all the glue you need.
#
# ```python
# run(["phenix.phaser", "hklin=obs.mtz", "labin=F=FOBS SIGF=SIGFOBS",
#      "model=6Y2E_search.pdb", "output.directory=phaser_out"])
# ```

# %% [markdown]
# ## 4. Evaluating a solution
#
# This is the section that actually matters, and all of it runs here.
#
# Phaser reports a translation-function Z score and a log-likelihood gain. Those numbers are
# useful but they are also opaque to a first-year student. The quantities underneath them are
# not opaque at all. A correct solution places atoms where the electron density is, so the
# calculated amplitudes track the observed amplitudes. A wrong solution does not.
#
# We build a correct solution by superposing the trimmed search model onto the deposited
# target, then measure how quickly the agreement decays as we push the model away from that
# position.

# %%
# Superpose the search model onto the target. This stands in for what Phaser found.
polymer_search = search[0][0].get_polymer()
polymer_target = target[0][0].get_polymer()

sup = gemmi.calculate_superposition(
    polymer_target, polymer_search, gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP
)
print(f"Superposition over {sup.count} CA pairs, RMSD {sup.rmsd:.2f} angstroms")

placed = search.clone()
placed[0].transform_pos_and_adp(sup.transform)
placed.cell = target.cell
placed.spacegroup_hm = target.spacegroup_hm
placed.setup_entities()
placed.write_pdb(str(WORK / "6Y2E_placed.pdb"))
print("Placed model written with the target cell and space group.")

# %%
D_MIN = 2.5  # Coarse enough to keep every cell in this section under a few seconds.


def as_dataset(miller_index, cell, spacegroup, **columns):
    """Build an rs.DataSet from a gemmi Miller array.

    The Miller columns have to carry the MTZ `HKL` dtype before they become the index,
    otherwise the DataSet cannot be written back out and will not join cleanly onto a
    DataSet that was read from a file. This is the one piece of bookkeeping that trips
    people up when they build a DataSet by hand.
    """
    idx = np.asarray(miller_index)
    ds_out = rs.DataSet({"H": idx[:, 0], "K": idx[:, 1], "L": idx[:, 2], **columns})
    for label in ("H", "K", "L"):
        ds_out[label] = ds_out[label].astype("HKL")
    ds_out = ds_out.set_index(["H", "K", "L"])
    ds_out.cell = cell
    ds_out.spacegroup = spacegroup
    ds_out.merged = True
    return ds_out


def calc_amplitudes(st, d_min=D_MIN):
    """Return a data frame of |Fcalc| indexed by Miller index, via an FFT of the model density."""
    dc = gemmi.DensityCalculatorX()
    dc.d_min = d_min
    dc.rate = 1.5
    dc.set_grid_cell_and_spacegroup(st)
    dc.put_model_density_on_grid(st[0])
    sf = gemmi.transform_map_to_f_phi(dc.grid, half_l=True)
    asu = sf.prepare_asu_data(dmin=d_min)
    return as_dataset(
        asu.miller_array,
        st.cell,
        gemmi.SpaceGroup(st.spacegroup_hm),
        FCALC=np.abs(asu.value_array),
    )


def r_factor_and_cc(st, observed=ds, d_min=D_MIN):
    """Scale |Fcalc| onto |Fobs| with one k and one B, then report R and the correlation."""
    calc = calc_amplitudes(st, d_min)
    merged = observed.join(calc, how="inner").dropna()
    if len(merged) < 50:
        return dict(n=len(merged), r=np.nan, cc=np.nan)

    fo = merged["FOBS"].to_numpy(dtype=float)
    fc = merged["FCALC"].to_numpy(dtype=float)
    merged = merged.compute_dHKL()
    s2 = 1.0 / merged["dHKL"].to_numpy(dtype=float) ** 2

    def residual(p):
        k, b = p
        return fo - k * np.exp(-b * s2 / 4.0) * fc

    k, b = scipy.optimize.least_squares(residual, [fo.mean() / max(fc.mean(), 1e-9), 0.0]).x
    fc_scaled = k * np.exp(-b * s2 / 4.0) * fc

    return dict(
        n=len(merged),
        k=float(k),
        B=float(b),
        r=float(np.abs(fo - fc_scaled).sum() / fo.sum()),
        cc=float(np.corrcoef(fo, fc_scaled)[0, 1]),
    )


baseline = r_factor_and_cc(placed)
print("Correct placement")
for key, value in baseline.items():
    print(f"   {key:4s} {value:.4f}" if isinstance(value, float) else f"   {key:4s} {value}")

# %% [markdown]
# An R factor near 0.35 to 0.45 and a correlation above 0.7 is what a correct, unrefined
# molecular replacement solution looks like. Anything near 0.55 is what random noise looks
# like. Give the class those two anchors before showing the next plot.

# %%
def perturb(st, axis, angle_deg, shift):
    """Rotate a structure about its own centroid and then translate it."""
    out = st.clone()
    pos = np.array([[a.pos.x, a.pos.y, a.pos.z] for m in out for c in m for r in c for a in r])
    centroid = pos.mean(axis=0)

    axis = np.asarray(axis, dtype=float)
    axis = axis / np.linalg.norm(axis)
    theta = np.deg2rad(angle_deg)
    kx, ky, kz = axis
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)

    for model in out:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    v = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                    w = R @ (v - centroid) + centroid + np.asarray(shift, dtype=float)
                    atom.pos = gemmi.Position(*w)
    return out


# How fast does the signal die as the model rotates away from the answer?
angles = [0, 1, 2, 3, 5, 8, 12, 20, 30, 45, 90]
rot_scan = pd.DataFrame(
    [dict(angle=a, **r_factor_and_cc(perturb(placed, [0, 1, 0], a, [0, 0, 0]))) for a in angles]
)
display(rot_scan[["angle", "n", "r", "cc"]])

# %%
# And as it slides away along a cell axis?
shifts = [0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0]
trans_scan = pd.DataFrame(
    [dict(shift=s, **r_factor_and_cc(perturb(placed, [0, 1, 0], 0, [s, 0, 0]))) for s in shifts]
)
display(trans_scan[["shift", "n", "r", "cc"]])

# %%
fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))

axes[0].plot(rot_scan["angle"], rot_scan["r"], marker="o", label="R factor")
axes[0].plot(rot_scan["angle"], rot_scan["cc"], marker="s", label="correlation")
axes[0].set_xlabel("rotation away from the solution, degrees")
axes[0].set_title("The rotation function has a narrow peak")

axes[1].plot(trans_scan["shift"], trans_scan["r"], marker="o", label="R factor")
axes[1].plot(trans_scan["shift"], trans_scan["cc"], marker="s", label="correlation")
axes[1].set_xlabel("translation away from the solution, angstroms")
axes[1].set_title("So does the translation function")

for ax in axes:
    ax.axhline(0.59, ls=":", c="grey")
    ax.annotate("random model", (ax.get_xlim()[1] * 0.45, 0.60), color="grey", fontsize=9)
    ax.set_ylabel("value")
    ax.legend()
    ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# %% [markdown]
# **The lesson.** Both curves are sharp. A model three degrees off, or two angstroms off, is
# already most of the way to worthless. That sharpness is why molecular replacement needs a
# systematic six-dimensional search rather than a hopeful guess, and it is why Phaser
# separates the rotation search from the translation search instead of searching all six
# dimensions at once.
#
# Ask the class to predict the width of the peak before you run the scan. Most will guess too
# wide.

# %% [markdown]
# ## 5. Initial refinement, one rigid body
#
# The first thing anyone does with a fresh MR solution is rigid-body refinement. Six
# parameters, no restraints, no geometry. Small enough to write from scratch, which makes it
# a good bridge into Demo 4 where a real refinement program takes over.

# %%
def rigid_body_refine(st, start_offset=(2.0, 0.0, 0.0), start_angle=3.0, maxiter=40):
    """Refine three rotations and three translations against the R factor."""
    history = []

    def objective(p):
        candidate = perturb(st, [1, 0, 0], p[0], [0, 0, 0])
        candidate = perturb(candidate, [0, 1, 0], p[1], [0, 0, 0])
        candidate = perturb(candidate, [0, 0, 1], p[2], p[3:6])
        r = r_factor_and_cc(candidate)["r"]
        history.append(dict(cycle=len(history), r=r, params=p.copy()))
        print(f"  cycle {len(history):3d}   R = {r:.4f}   {np.round(p, 3)}")
        return r

    start = np.array([start_angle, 0.0, 0.0, *start_offset])
    print("Starting from a deliberately displaced model.")
    print(f"  start        R = {objective(start):.4f}")
    result = scipy.optimize.minimize(
        objective, start, method="Powell", options=dict(maxiter=maxiter, xtol=0.05, ftol=1e-4)
    )
    return result, pd.DataFrame(history)


result, history = rigid_body_refine(placed)
print("\nRefined parameters :", np.round(result.x, 3))
print("Final R factor     :", round(result.fun, 4))
print("Target R factor    :", round(baseline["r"], 4))

# %%
fig, ax = plt.subplots(figsize=(7, 4.2))
ax.plot(history["cycle"], history["r"], marker=".", lw=1)
ax.axhline(baseline["r"], ls="--", c="tab:green", label="correctly placed model")
ax.set_xlabel("objective evaluation")
ax.set_ylabel("R factor")
ax.set_title("Rigid-body refinement pulls a displaced model back onto the answer")
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# %% [markdown]
# ## 6. What the real programs add
#
# The hand-rolled refinement above optimizes one number. A production molecular replacement
# program does considerably more, and naming the differences keeps students from mistaking
# the demonstration for the real thing.
#
# * **A likelihood target.** Phaser maximizes a log-likelihood that models the error in the
#   search model explicitly, rather than minimizing an R factor that treats every reflection
#   as equally trustworthy.
# * **A systematic search.** The fast rotation function scores thousands of orientations
#   before any translation is attempted.
# * **Packing checks.** A solution that puts two copies of the protein in the same place is
#   rejected before it is reported.
# * **Space group ambiguity.** `sgalternative=all` tries every space group consistent with
#   the observed systematic absences, which routinely resolves an enantiomorph the data alone
#   cannot distinguish.
# * **Multiple copies.** Real crystals often have several molecules per asymmetric unit, found
#   one at a time with each previous copy held fixed.

# %% [markdown]
# ## Exercises
#
# 1. Rerun the rotation scan with the poly-alanine model. Does the peak get wider or narrower?
#    Explain the result in terms of how much of the scattering the model accounts for.
# 2. Repeat the scan at 3.5 angstroms and at 2.0 angstroms. How does resolution change the
#    width of the peak?
# 3. Compute the Matthews coefficient for this cell and space group, assuming 306 residues per
#    chain, and decide how many copies belong in the asymmetric unit. Compare against 6YNQ.
# 4. Replace the R factor in `rigid_body_refine` with a correlation on intensities and see
#    whether the refinement converges from a larger starting displacement.
# 5. Run the real PHASER on a workstation with `6Y2E_polyala.pdb` and compare the reported
#    TFZ against the R factor you get here for the same placement.
#
# ## What to carry forward
#
# `6Y2E_placed.pdb` and `obs.mtz` are the inputs to Demo 4.
#
# ## Sources
#
# * CHAINSAW. Stein, N. (2008). *J. Appl. Cryst.* **41**, 641-643. https://doi.org/10.1107/S0021889808006985
# * Phaser. McCoy, A. J. et al. (2007). *J. Appl. Cryst.* **40**, 658-674. https://doi.org/10.1107/S0021889807021206
# * MOLREP. Vagin, A. and Teplyakov, A. (2010). *Acta Cryst.* **D66**, 22-25. https://doi.org/10.1107/S0907444909042589
# * Gemmi. https://gemmi.readthedocs.io
# * PDB entries 6YNQ and 6Y2E. https://www.rcsb.org
