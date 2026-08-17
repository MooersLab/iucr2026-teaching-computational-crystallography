# Live demonstration notebooks

Companion notebooks for **Teaching Crystallographic Computing with Computational Notebooks**,
Teaching Computational Crystallography Minisymposium, 27th Congress and General Assembly of
the IUCr, Calgary, Alberta, 2026 August 17.

Blaine Mooers, PhD
Department of Biochemistry and Physiology, University of Oklahoma Health Campus
blaine-mooers@ou.edu

## The six demonstrations

`main.pdf` has six slides that read **LIVE DEMONSTRATION** and say *Switching to Jupyter
notebook workspace*. One notebook covers each.

| Slide | Frame title | Notebook | Runs in Colab |
|------:|-------------|----------|---------------|
| 7 | DIALS in Jupyter | `demo1_DIALS_in_Jupyter.ipynb` | Yes, after a conda install |
| 8 | XDS via Subprocess Calls | `demo2_XDS_via_subprocess.ipynb` | Yes, with your XDS license link |
| 10 | Molecular Replacement Workflow | `demo3_molecular_replacement.ipynb` | Partly. See below |
| 12 | Refinement Workflow | `demo4_refinement_workflow.ipynb` | Partly. See below |
| 14 | reciprocalspaceship and Gemmi | `demo5_reciprocalspaceship_gemmi.ipynb` | Yes, entirely |
| 18 | Reducing Errors with Snippets | `demo6_snippet_libraries.ipynb` | Yes, entirely |

Every notebook lists the same four bullet points its slide lists, in the same order, so you
can move from the slide to the notebook without losing the thread.

## What "partly" means

Demos 3 and 4 name software that cannot be installed inside a free Colab session. CCP4 is six
to eight gigabytes with a graphical front end. Phenix ships as a personalized installer and is
on no public package index.

Rather than show dead cells, those two notebooks split the work.

* The MOLREP, PHASER, CHAINSAW, `phenix.refine`, and REFMAC5 command lines appear exactly as
  you would type them on a workstation.
* Everything a student needs to *understand* those programs runs live on `gemmi`, `numpy`,
  and `scipy`. Demo 3 measures how fast the molecular replacement signal decays as a model
  moves away from the answer. Demo 4 computes R-work and R-free from scratch and plots the
  gap between them widening.

The parts that run are the parts worth class time. Launching Phaser is one command. Knowing
whether the answer is right is the lesson.

## Data

Every notebook uses the same crystal, so the six demonstrations tell one story.

* **Raw frames.** Zenodo record [3752540](https://doi.org/10.5281/zenodo.3752540). 686 CBF
  frames of SARS-CoV-2 main protease bound to 2-methyl-1-tetralone, collected on beamline P11
  at PETRA III. 0.2 degrees per frame, 40 ms exposure, wavelength 1.033214 angstroms, detector
  at 200 mm. Each frame is 6.2 MB, so the full sweep is 6.2 GB. The notebooks take a
  60-frame wedge, which is twelve degrees and about 370 MB.
* **Target structure.** PDB entry [6YNQ](https://www.rcsb.org/structure/6YNQ).
* **Search model.** PDB entry [6Y2E](https://www.rcsb.org/structure/6Y2E), the apo enzyme.

## How the notebooks chain together

```
Demo 1  frames ──> scaled.mtz ──┐
Demo 2  frames ──> xds.mtz   ───┤ compared against each other
                                │
Demo 3  scaled.mtz + 6Y2E ──> 6Y2E_placed.pdb
                                │
Demo 4  6Y2E_placed.pdb + scaled.mtz ──> R-work, R-free, geometry

Demo 5  stands alone. Deposited data for 6YNQ.
Demo 6  stands alone. No data at all.
```

Each notebook also runs on its own. When an upstream output is missing, the notebook falls
back to the deposited data for 6YNQ and says so in its output.

## Running order in a talk

If you have time for only one, show **Demo 5**. It installs in fifteen seconds, never fails,
and makes the case that a student who knows pandas already knows most of the library.

If you have time for two, add **Demo 6**. It needs no data and no network, and the d-spacing
comparison in section 5 is the argument that changes minds.

Demos 1 and 2 need to be prebaked. Run them the night before, save a copy to Google Drive,
and open that copy in class. Cells marked **PREBAKE** are the slow ones.

## Files

Each notebook is paired with a Jupytext percent-format `.py` file of the same name. The
notebook is what you open. The `.py` file is what you commit, because it diffs cleanly.

To regenerate a notebook after editing the `.py`:

```bash
jupytext --to ipynb demo1_DIALS_in_Jupyter.py
```

To keep the pair in sync automatically inside Jupyter, install `jupytext` and the pairing in
each notebook's metadata will do the rest.

```bash
pip install jupytext
jupytext --set-formats ipynb,py:percent demo1_DIALS_in_Jupyter.ipynb
```

## Dependencies

Demos 5 and 6 need only this.

```bash
pip install gemmi reciprocalspaceship pandas numpy scipy matplotlib
```

Demos 3 and 4 need the same list. Demo 1 additionally needs DIALS, which comes from
conda-forge.

```bash
conda create -n dials -c conda-forge dials gemmi reciprocalspaceship matplotlib pandas
```

Demo 2 needs an XDS binary, which requires an academic license from
[xds.mr.mpg.de](https://xds.mr.mpg.de). Paste your download link into the `XDS_URL` variable
in section 2.

## Exercises

Every notebook closes with four or five exercises. They are written to be assigned, not just
read. Several ask the student to break something on purpose and explain the damage, which is
the fastest way to build the instinct that a plausible number is not a correct one.
