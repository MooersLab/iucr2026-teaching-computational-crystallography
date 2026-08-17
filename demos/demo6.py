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
# # Demo 6. Reducing errors with code snippet libraries
#
# **Companion to slide 18 of *Teaching Crystallographic Computing with Computational Notebooks*.**
# 27th Congress and General Assembly of the IUCr, Calgary, Alberta, 2026 August 17.
#
# Blaine Mooers, PhD. Department of Biochemistry and Physiology, University of Oklahoma Health Campus.
# blaine-mooers@ou.edu
#
# ## What this demonstration covers
#
# 1. Code snippet libraries as software engineering tools.
# 2. GhostText for tab triggers in notebooks.
# 3. Advantages over AI-generated code for beginners.
# 4. Building personal snippet collections.
#
# ## The argument in one paragraph
#
# A first-year graduate student who has never programmed faces two problems at once. The
# first is crystallography. The second is Python. A snippet library removes the second problem
# for the twenty operations that recur in every session, so the student spends their attention
# on the science. The snippets are tested once, by you, and then reused without variation.
# That property is the whole value, and it is the property a language model cannot offer.

# %%
import os

os.environ["PYDEVD_DISABLE_FILE_VALIDATION"] = "1"

import pathlib
import sys

IN_COLAB = "google.colab" in sys.modules or os.path.isdir("/content")
WORK = pathlib.Path("/content/snippets" if IN_COLAB else "./snippets").resolve()
WORK.mkdir(parents=True, exist_ok=True)
os.chdir(WORK)
print("Working in", WORK)

# %%
# !pip -q install gemmi numpy pandas

# %%
import json
import re
import textwrap

import gemmi
import numpy as np
import pandas as pd

print("gemmi :", gemmi.__version__)

# %% [markdown]
# ## 1. What a snippet actually is
#
# A snippet is a named block of code with numbered tab stops. The syntax below is the
# LSP-standard form that VS Code, YASnippet, UltiSnips, and every modern editor understand.
#
# * `${1:default}` is the first tab stop, with a default value.
# * `${2:default}` is the second, and so on.
# * `$0` is where the cursor lands when you are done.
# * The same number used twice is a mirror, updated as you type.
#
# The whole idea fits in a dictionary.

# %%
SNIPPETS = {
    "rsread": {
        "prefix": "rsread",
        "description": "Read an MTZ into a reciprocalspaceship DataSet and report its metadata",
        "requires": "",
        "produces": "ds (a reciprocalspaceship DataSet with a dHKL column)",
        "body": textwrap.dedent(
            """\
            import reciprocalspaceship as rs

            ${1:ds} = rs.read_mtz("${2:data.mtz}")
            ${1:ds} = ${1:ds}.compute_dHKL()
            print("Cell        :", ${1:ds}.cell)
            print("Space group :", ${1:ds}.spacegroup.hm)
            print("Resolution  : {:.2f} to {:.2f} angstroms".format(
                ${1:ds}["dHKL"].max(), ${1:ds}["dHKL"].min()))
            ${1:ds}.head()$0"""
        ),
    },
    "gemread": {
        "prefix": "gemread",
        "description": "Read a coordinate file with gemmi and clean it up",
        "requires": "",
        "produces": "st (a gemmi Structure with waters and hydrogens removed)",
        "body": textwrap.dedent(
            """\
            import gemmi

            ${1:st} = gemmi.read_structure("${2:model.cif}")
            ${1:st}.setup_entities()
            ${1:st}.remove_alternative_conformations()
            ${1:st}.remove_hydrogens()
            print("${1:st}:", ${1:st}[0].count_atom_sites(), "atoms in",
                  len(${1:st}[0]), "chains")$0"""
        ),
    },
    "dspacing": {
        "prefix": "dspacing",
        "description": "Resolution of a Miller index for any lattice, via gemmi",
        "requires": "",
        "produces": "d (the d-spacing in angstroms for the given cell and hkl)",
        "body": textwrap.dedent(
            """\
            import gemmi

            cell = gemmi.UnitCell(${1:a}, ${2:b}, ${3:c}, ${4:alpha}, ${5:beta}, ${6:gamma})
            d = cell.calculate_d((${7:h}, ${8:k}, ${9:l}))
            print(f"d = {d:.3f} angstroms")$0"""
        ),
    },
    "resbins": {
        "prefix": "resbins",
        "description": "Equal-population resolution shells with a groupby summary",
        "requires": "ds (a reciprocalspaceship DataSet with a dHKL column, from rsread)",
        "produces": "summary (a DataFrame keyed by shell with d_low, d_high, n)",
        "body": textwrap.dedent(
            """\
            import numpy as np

            ${1:ds}["invd2"] = 1.0 / ${1:ds}["dHKL"] ** 2
            edges = np.quantile(${1:ds}["invd2"], np.linspace(0, 1, ${2:21}))
            ${1:ds}["shell"] = np.clip(np.digitize(${1:ds}["invd2"], edges[1:-1]), 0, ${2:21} - 2)
            summary = ${1:ds}.groupby("shell").agg(
                d_low=("dHKL", "max"), d_high=("dHKL", "min"), n=("dHKL", "size"))
            summary$0"""
        ),
    },
    "runext": {
        "prefix": "runext",
        "description": "Run an external crystallographic program and capture everything",
        "requires": "",
        "produces": "proc (a subprocess.CompletedProcess with stdout, stderr, returncode)",
        "body": textwrap.dedent(
            """\
            import subprocess

            proc = subprocess.run(
                ${1:["dials.version"]},
                capture_output=True, text=True, timeout=${2:3600})
            print("return code", proc.returncode)
            print(proc.stdout[-2000:])
            if proc.returncode != 0:
                raise RuntimeError(proc.stderr[-2000:])$0"""
        ),
    },
}

print(f"{len(SNIPPETS)} snippets in the library\n")
for name, snip in SNIPPETS.items():
    print(f"  {name:10s} {snip['description']}")

# %% [markdown]
# ## 2. Expanding a snippet
#
# In an editor, pressing Tab does this. In a notebook without an editor extension, forty
# lines of Python do the same thing, and writing those forty lines in front of a class
# removes the mystery from the editor feature.

# %%
TABSTOP = re.compile(r"\$\{(\d+):([^}]*)\}|\$(\d+)")


def tab_stops(body):
    """Return the ordered tab stops of a snippet with their default values."""
    found = {}
    for m in TABSTOP.finditer(body):
        number = int(m.group(1) or m.group(3))
        default = m.group(2) if m.group(1) else ""
        found.setdefault(number, default)
    return dict(sorted(found.items()))


def defaults_summary(body):
    """Format the default values of a snippet's tab stops for menu display.

    Returns a string like '$1=ds, $2=data.mtz'. Skips $0 (the final cursor
    position). Bare $1-style stops with no default show as '$1' alone.
    """
    parts = []
    for n, default in tab_stops(body).items():
        if n == 0:
            continue
        parts.append(f"${n}={default}" if default else f"${n}")
    return ", ".join(parts)


def expand(name, *values, **named):
    """Fill a snippet's tab stops, positionally or by number."""
    body = SNIPPETS[name]["body"]
    stops = tab_stops(body)
    filled = dict(stops)
    for i, value in enumerate(values, start=1):
        filled[i] = str(value)
    for key, value in named.items():
        filled[int(key.lstrip("s"))] = str(value)

    def substitute(m):
        number = int(m.group(1) or m.group(3))
        return "" if number == 0 else filled.get(number, "")

    return TABSTOP.sub(substitute, body)


print("Tab stops in 'dspacing':", tab_stops(SNIPPETS["dspacing"]["body"]))
print()
print(expand("dspacing", 45.11, 53.98, 114.79, 90, 101.5, 90, 2, 1, 0))

# %% [markdown]
# ## 3. Tab triggers inside a notebook, with no extension
#
# IPython lets a magic write into the next cell. That is enough to build a working tab
# trigger, and it works in Colab, in JupyterLab, and in the classic notebook without
# installing anything.

# %%
from IPython.core.magic import register_line_magic


@register_line_magic
def snip(line):
    """Insert a snippet into a new cell below. Usage: %snip rsread 6ynq.mtz"""
    parts = line.split()
    if not parts:
        for name, s in SNIPPETS.items():
            print(f"  %snip {name:10s} {s['description']}")
            print(f"  {'':10s}   defaults: {defaults_summary(s['body'])}")
            if s.get("requires"):
                print(f"  {'':10s}   requires: {s['requires']}")
            if s.get("produces"):
                print(f"  {'':10s}   produces: {s['produces']}")
        return
    name, args = parts[0], parts[1:]
    if name not in SNIPPETS:
        print(f"No snippet named {name}. Run %snip with no arguments for the list.")
        return
    code = expand(name, *args)
    ip = get_ipython()  # noqa: F821
    ip.set_next_input(f"# snippet: {name}\n{code}", replace=False)
    print(f"Inserted '{name}' into a new cell below.")


print("Magic registered. Try:  %snip")

# %%
# Run this cell, then look at the cell that appears below it.
# %snip rsread 6ynq_data 6YNQ.mtz

# %%
# The library, browsable from the notebook.
display(
    pd.DataFrame(
        [
            dict(trigger=name,
                 description=s["description"],
                 requires=s.get("requires", ""),
                 produces=s.get("produces", ""),
                 defaults=defaults_summary(s["body"]),
                 stops=len(tab_stops(s["body"])),
                 lines=s["body"].count("\n") + 1)
            for name, s in SNIPPETS.items()
        ]
    )
)

# %% [markdown]
# ## 4. GhostText, for the students who already have an editor
#
# The `%snip` magic above is a teaching device. The production answer is GhostText, which
# connects a text area in the browser to a real editor. You edit a Jupyter cell in Emacs, Vim,
# VS Code, or Sublime, and every keystroke appears in the browser. Every snippet system,
# keybinding, and linter you already configured comes with you.
#
# The setup is three steps.
#
# 1. Install the GhostText browser extension in Chrome or Firefox.
# 2. Install the matching editor package. `atomic-chrome` for Emacs, `ghost-text.vim` for Vim,
#    `GhostText` from the marketplace for VS Code.
# 3. Start the editor server. In Emacs that is `M-x atomic-chrome-start-server`.
#
# Then click into a notebook cell and press the GhostText button. Tab triggers, tab stops, and
# mirrored fields all work, because the editing is happening in the editor.
#
# **Why this matters more than it sounds.** JupyterLab extensions break on every major
# release, and an instructor who depends on one will eventually stand in front of a class with
# a broken environment. GhostText moves the fragile part out of the notebook stack entirely.
# The browser extension and the editor package are versioned independently of Jupyter.
#
# Snippet libraries for structural biology, in editor-native formats, live at
# [github.com/MooersLab](https://github.com/MooersLab). The `cctbxsnips` and `pymolsnips`
# collections cover the two APIs that come up most in a crystallography class.

# %% [markdown]
# ## 5. Why a tested snippet beats generated code for a beginner
#
# This is the part of the demonstration that changes minds, because it is a measurement rather
# than an opinion.
#
# The formula below is what a plausible answer to "how do I compute the resolution of a
# reflection?" looks like. It is correct. It is also correct only for orthorhombic lattices,
# and it fails silently everywhere else. A student who cannot yet read the formula has no way
# to detect the failure, and the wrong numbers propagate into a resolution cutoff, a shell
# table, and a figure.

# %%
def d_spacing_naive(a, b, c, h, k, l):
    """The formula a beginner is most likely to be handed. Orthorhombic only."""
    return 1.0 / np.sqrt((h / a) ** 2 + (k / b) ** 2 + (l / c) ** 2)


# The Mpro crystal is monoclinic C2. Beta is not 90 degrees.
CELL = gemmi.UnitCell(45.11, 53.98, 114.79, 90.0, 101.5, 90.0)
print("Test cell :", CELL, "\n")

reflections = [(2, 0, 0), (0, 2, 0), (0, 0, 4), (1, 1, 1), (4, 2, 6), (8, 0, 12), (2, 4, 20)]

rows = []
for h, k, l in reflections:
    truth = CELL.calculate_d((h, k, l))
    naive = d_spacing_naive(CELL.a, CELL.b, CELL.c, h, k, l)
    rows.append(
        dict(hkl=f"{h} {k} {l}", correct=truth, naive=naive,
             error_pct=100 * (naive - truth) / truth)
    )

comparison = pd.DataFrame(rows)
display(comparison.round(3))
print(f"\nWorst error : {comparison['error_pct'].abs().max():.1f} percent")
print("Every number in the 'naive' column looks reasonable. That is the problem.")

# %% [markdown]
# ### The same comparison, as a property
#
# A snippet has a property that generated code does not. Expand it a hundred times and you get
# the same hundred characters. That is not a claim about quality, it is a claim about variance,
# and variance is what makes a beginner's environment unteachable.

# %%
expansions = {expand("dspacing", 45.11, 53.98, 114.79, 90, 101.5, 90, 2, 1, 0) for _ in range(100)}
noun = "result" if len(expansions) == 1 else "results"
print(f"100 expansions produced {len(expansions)} distinct {noun}.")

# %% [markdown]
# **Where each tool belongs.** This is not an argument against language models. It is an
# argument about sequencing.
#
# | | Snippet library | Generated code |
# |---|---|---|
# | Output for the same request | Identical every time | Varies |
# | Correctness | Established once, by you | Must be checked every time |
# | Who can check it | Anyone, because you already did | Only someone who already knows the answer |
# | Covers a novel problem | No | Yes |
# | Teaches the API | Yes, by repeated exposure | Weakly, because the code arrives finished |
# | Works offline, in a locked-down lab | Yes | Usually not |
#
# Give beginners snippets for the recurring twenty operations. Bring in generation once they
# can read the output well enough to reject it. A student who cannot yet tell a correct
# d-spacing from an incorrect one is not in a position to supervise a model, and supervision
# is the entire job.

# %% [markdown]
# ## 6. Building and exporting a personal collection
#
# A collection grows one snippet at a time, out of the code you just finished debugging. The
# rule that makes it work is to add the snippet at the moment it works, not later.

# %%
def add_snippet(name, description, body, store=SNIPPETS):
    """Add a snippet after checking that its tab stops are numbered sensibly."""
    stops = tab_stops(body)
    numbers = [n for n in stops if n != 0]
    if numbers and sorted(numbers) != list(range(1, max(numbers) + 1)):
        raise ValueError(f"Tab stops must run 1..n with no gaps. Got {sorted(numbers)}.")
    store[name] = dict(prefix=name, description=description, body=body)
    print(f"Added '{name}' with {len(numbers)} tab stops.")
    return store[name]


add_snippet(
    "wilsonb",
    "Wilson B factor from a DataSet that already has dHKL",
    textwrap.dedent(
        """\
        import numpy as np

        w = ${1:ds}[["${2:I}", "dHKL"]].dropna().copy()
        w["invd2"] = 1.0 / w["dHKL"] ** 2
        edges = np.quantile(w["invd2"], np.linspace(0, 1, 21))
        w["shell"] = np.clip(np.digitize(w["invd2"], edges[1:-1]), 0, 19)
        shells = w.groupby("shell").agg(invd2=("invd2", "mean"),
                                        d=("dHKL", "mean"),
                                        Imean=("${2:I}", "mean"))
        high = shells[shells["d"] < ${3:4.0}]
        slope, _ = np.polyfit(high["invd2"], np.log(high["Imean"]), 1)
        print(f"Wilson B = {-2 * slope:.1f} square angstroms")$0"""
    ),
)

# %%
def to_vscode(store=SNIPPETS):
    """Export to the VS Code snippet JSON format."""
    return json.dumps(
        {
            s["description"]: {"prefix": s["prefix"], "body": s["body"].split("\n"),
                               "description": s["description"]}
            for s in store.values()
        },
        indent=2,
    )


def to_ultisnips(store=SNIPPETS):
    """Export to the UltiSnips format used by Vim and Neovim."""
    blocks = []
    for s in store.values():
        blocks.append(f'snippet {s["prefix"]} "{s["description"]}" b\n{s["body"]}\nendsnippet\n')
    return "\n".join(blocks)


def to_yasnippet(name, store=SNIPPETS):
    """Export one snippet to the YASnippet format used by Emacs. One file per snippet."""
    s = store[name]
    return f"# -*- mode: snippet -*-\n# name: {s['description']}\n# key: {s['prefix']}\n# --\n{s['body']}\n"


pathlib.Path("crystallography.code-snippets").write_text(to_vscode())
pathlib.Path("crystallography.snippets").write_text(to_ultisnips())
yas_dir = pathlib.Path("yasnippets/python-mode")
yas_dir.mkdir(parents=True, exist_ok=True)
for name in SNIPPETS:
    (yas_dir / name).write_text(to_yasnippet(name))

print("Exported to three editor formats:")
for path in sorted(pathlib.Path(".").rglob("*")):
    if path.is_file():
        print(f"   {path}  ({path.stat().st_size} bytes)")

# %%
print(to_vscode()[:900], "\n...")

# %% [markdown]
# ## 7. Handing the collection to students
#
# One repository, one install line per editor, and a table of triggers on the first page.
#
# * **VS Code.** Drop the `.code-snippets` file into `.vscode/` in the project folder.
# * **Vim or Neovim.** Put the `.snippets` file where UltiSnips looks, usually
#   `~/.vim/UltiSnips/python.snippets`.
# * **Emacs.** Put the `yasnippets/python-mode/` folder on `yas-snippet-dirs`.
# * **JupyterLab or Colab.** Use GhostText into any of the above, or ship the `%snip` magic in
#   a small module the students import at the top of every notebook.
#
# Ask each student to add one snippet per week from their own debugged code. By the end of a
# semester the class has built a shared library, and building it taught them more than
# receiving it would have.

# %% [markdown]
# ## Exercises
#
# 1. Add a snippet for the DIALS import and spot-finding pair from Demo 1. Give it three tab
#    stops.
# 2. Extend `expand` to support mirrored tab stops that appear more than twice, then verify it
#    on the `rsread` snippet.
# 3. Write a test that expands every snippet with its default values and byte-compiles the
#    result. Run it on the whole library. This is the check that keeps a collection honest.
# 4. Take a piece of crystallographic code a language model wrote for you last week and find
#    the assumption it made silently. Turn the corrected version into a snippet.
# 5. Add an export function for your own editor if it is not one of the three above.
#
# ## Sources
#
# * GhostText. https://ghosttext.fregante.com
# * MooersLab snippet libraries. https://github.com/MooersLab
# * Gemmi. https://gemmi.readthedocs.io
