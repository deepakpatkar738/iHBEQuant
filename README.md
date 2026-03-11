<p align="center">
  <img src="iHBEQuant_logo.png" alt="iHBEQuant Logo" width="320"/>
</p>

<h1 align="center">iHBEQuant — Individual Hydrogen Bond Energy Quantifier</h1>

<p align="center">
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.7%2B-blue.svg" alt="Python 3.7+"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
  <img src="https://img.shields.io/badge/version-1.2-orange.svg" alt="Version">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Windows%20%7C%20macOS-lightgrey.svg" alt="Platform">
</p>

**iHBEQuant** is a Python tool that automates the quantification of individual hydrogen bond (HB) energies in molecular clusters using the **Molecular Tailoring Approach (MTA)**. Given a cluster geometry in XYZ format, it automatically detects hydrogen bonds, generates all required Gaussian input files, optionally submits and monitors calculations, and parses the results into a clean energy summary in a plain-text file.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output](#output)
- [Methodology](#methodology)
- [Example](#example)
- [Citations](#citations)
- [Author](#author)
- [License](#license)

---

## Overview

Hydrogen bonds are central to molecular recognition, crystal packing, protein folding, and the bulk properties of liquids and solids. Standard quantum chemical calculations yield only a *total* interaction energy for a cluster, making it impossible to directly assess the role of each individual hydrogen bond.

**iHBEQuant** solves this by implementing the **Molecular Tailoring Approach (MTA)** for individual hydrogen bond energy (HBE) quantification. In MTA, the cluster is systematically fragmented — by removing the donor monomer, the acceptor monomer, or both — and the individual HBE is recovered from the total and fragment energies without any empirical parameters. iHBEQuant fully automates this fragmentation, Gaussian job management, energy parsing, and reporting workflow.

---

## Features

- **Automatic HB detection** — geometry-driven D–H···A search with fully configurable distance and angle thresholds
- **Clear file-naming convention** — all files written into a per-system subdirectory:
  - Parent molecule: `F1O2_M.xyz`, `F1O2_M.gjf`
  - Monomers: `F1O2_mono_M1.gjf`, `_mono_M2.gjf`, …
  - Primary fragments: `HB1_F1O2_frag_M1.gjf` (donor side), `HB1_F1O2_frag_M3.gjf` (acceptor side)
  - Overlap fragment: `HB1_F1O2_frag_M1M3.gjf`
  - Dimer: `HB1_F1O2_dimer_M1M3.gjf`
- **Three run modes** — `gen` (write files only), `run` (submit to Gaussian 16), `cal` (parse existing logs)
- **MTA, dimer, and cooperative energies** — decomposes each HBE into in-cluster MTA, pairwise dimer, and cooperativity contributions
- **Binding energy and MTA error** — total cluster binding energy cross-validated against the sum of individual HBEs
- **Coordinates embedded in summary** — all fragment geometries written directly into `_summary.txt`
- **Multi-file batch processing** — space-separated `.xyz` list in `INPUT.cfg`, processed sequentially
- **Supports DFT, MP2, CCSD(T)** — any single-point method available in Gaussian 16
- **Zero dependencies** — Python standard library only; no `pip install` required

---

## Requirements

- **Python 3.7+**
- **[Gaussian 16](https://gaussian.com/)** with `g16` on the system PATH — required only for `run` and `cal` modes; `gen` mode works without Gaussian

---

## Installation

```bash
git clone https://github.com/deepakpatkar738/iHBEQuant.git
cd iHBEQuant
```

No compilation or package installation needed:

```bash
python iHBEQuant.py
```

---

## Quick Start

1. Place your cluster geometry as a standard XYZ file (e.g., `F1O2.xyz`) in the working directory.
2. Edit `INPUT.cfg` — set the XYZ filename, method, basis set, charge, multiplicity, and run mode.
3. Run:

```bash
python iHBEQuant.py
```

Results appear in a subdirectory named after your XYZ file (e.g., `F1O2/`).

---

## Usage

```bash
python iHBEQuant.py                  # use default INPUT.cfg
python iHBEQuant.py -c myrun.cfg     # custom config file
python iHBEQuant.py -c myrun.cfg -v  # verbose: print all .gjf filenames
```

### Command-Line Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-c`, `--config` | `INPUT.cfg` | Path to the configuration file |
| `-v`, `--verbose` | off | Print all generated `.gjf` filenames to screen |

---

## Configuration

All run parameters are controlled through a plain-text `INPUT.cfg` file. Inline `#` comments are stripped automatically.

```ini
[SYSTEM]
xyz_file  = F1O2.xyz        # space-separated list of .xyz files
nproc     = 8               # CPU cores for Gaussian
mem       = 8GB             # memory for Gaussian

[THEORY]
method       = MP2          # DFT functional, MP2, CCSD(T), etc.
basis        = aug-cc-pVTZ  # basis set
keywords     = scf=tight    # extra Gaussian route keywords (appended to every .gjf)
charge       = 0            # total charge of the full cluster
multiplicity = 1            # spin multiplicity (1 = singlet, 2 = doublet, ...)

[THRESHOLDS]
hb_min_dist  = 1.2          # H···A minimum distance (Å)
hb_max_dist  = 2.8          # H···A maximum distance (Å)
hb_min_angle = 100.0        # D–H···A minimum angle (degrees)
hb_max_angle = 180.0        # D–H···A maximum angle (degrees)
neighbor_th  = 2.5          # covalent neighbor cutoff (Å)

[RUN]
mode = gen
# gen  →  write .gjf input files only
# run  →  write .gjf files + submit to g16
# cal  →  parse existing .log files only
```

> **Important:** `charge` and `multiplicity` apply to the full cluster and are inherited by all fragment calculations. iHBEQuant does **not** auto-detect these from the XYZ file.

---

## Output

For each input XYZ file (e.g., `F1O2.xyz`), a subdirectory `F1O2/` is created containing:

### Gaussian Input Files

| File | Role |
|------|------|
| `F1O2_M.xyz` | Copy of the input geometry (parent molecule reference) |
| `F1O2_M.gjf` | Full N-monomer cluster single-point |
| `F1O2_mono_M1.gjf`, `_mono_M2.gjf`, … | Isolated monomers |
| `HB1_F1O2_frag_M1.gjf` | Primary fragment — donor side (acceptor removed) |
| `HB1_F1O2_frag_M3.gjf` | Primary fragment — acceptor side (donor removed) |
| `HB1_F1O2_frag_M1M3.gjf` | Overlap fragment — both donor and acceptor removed |
| `HB1_F1O2_dimer_M1M3.gjf` | Isolated donor–acceptor dimer |

Files prefixed `HB1_` are replicated as `HB2_`, `HB3_`, … for each detected hydrogen bond.

### Summary File

`F1O2_summary.txt` contains the run header (method, thresholds, statistics), the HB energy table, per-HB calculation details with raw Hartree energies, and all fragment coordinates in XYZ format.

---

## Methodology

iHBEQuant implements the **Molecular Tailoring Approach** for individual hydrogen bond energies as formulated by Ahirwar, Gadre, and Deshmukh (2020) and extended by Patkar et al. (2021–2022).

For HB number $n$ between donor monomer **D** (index $d$) and acceptor monomer **A** (index $a$) in an $N$-monomer cluster, seven single-point energies are required:

### Fragment Definitions

| Symbol | File generated | Description |
|--------|---------------|-------------|
| $E_{\mathrm{M}}$ | F1O2\_M.gjf | Full $N$-monomer parent cluster |
| $E_{\mathrm{frag\text{-}D}}$ | HBn\_F1O2\_frag\_M*d*.gjf | Cluster with acceptor removed (donor side) |
| $E_{\mathrm{frag\text{-}A}}$ | HBn\_F1O2\_frag\_M*a*.gjf | Cluster with donor removed (acceptor side) |
| $E_{\mathrm{frag\text{-}DA}}$ | HBn\_F1O2\_frag\_M*d*_M*a*.gjf | Cluster with both donor and acceptor removed |
| $E_{\mathrm{dimer}}$ | HBn\_F1O2\_dimer\_M*d*_M*a*.gjf | Isolated donor–acceptor dimer |
| $E_{\mathrm{mono\text{-}D}}$ | F1O2\_mono\_M*d*.gjf | Isolated donor monomer |
| $E_{\mathrm{mono\text{-}A}}$ | F1O2\_mono\_M*a*.gjf | Isolated acceptor monomer |

---

### 1. MTA Hydrogen Bond Energy

The individual HB energy of HB $n$ within the full cluster environment:

$$\boxed{
E_{\mathrm{MTA\text{-}HB}n} \;=\;
\Bigl[\,\bigl(E_{\mathrm{frag\text{-}D}} + E_{\mathrm{frag\text{-}A}} - E_{\mathrm{frag\text{-}DA}}\bigr)
\;-\ E_{\mathrm{M}}\,\Bigr]
\;\times\ 627.51 \;\;\text{kcal mol}^{-1}
}$$

The three-fragment inclusion–exclusion $(E_{\mathrm{frag\text{-}D}} + E_{\mathrm{frag\text{-}A}} - E_{\mathrm{frag\text{-}DA}})$ reconstructs the energy contribution of the D···A pair within the cluster, and subtracting $E_{\mathrm{M}}$ isolates the individual HB contribution.

---

### 2. Dimer Hydrogen Bond Energy

The pairwise interaction energy of the isolated donor–acceptor pair in vacuum:

$$\boxed{
E_{\mathrm{Dimer\text{-}HB}n} \;=\;
\Bigl[\,E_{\mathrm{mono\text{-}D}} + E_{\mathrm{mono\text{-}A}} - E_{\mathrm{dimer}}\,\Bigr]
\;\times\; 627.51 \;\;\text{kcal mol}^{-1}
}$$

This is the conventional counterpoise-free interaction energy of the donor–acceptor pair with no influence from the remaining monomers.

---

### 3. Cooperativity Energy

The cooperative enhancement (positive) or anti-cooperative weakening (negative) arising from the surrounding cluster environment:

$$\boxed{
E_{\mathrm{Coop\text{-}HB}n} \;=\; E_{\mathrm{MTA\text{-}HB}n} \;-\; E_{\mathrm{Dimer\text{-}HB}n}
}$$

| Sign | Interpretation |
|------|---------------|
| $E_{\mathrm{Coop}} > 0$ | Cluster environment **strengthens** the HB (cooperative) |
| $E_{\mathrm{Coop}} < 0$ | Cluster environment **weakens** the HB (anti-cooperative) |
| $E_{\mathrm{Coop}} = 0$ | No environmental influence (purely pairwise) |

---

### 4. Total Cluster Binding Energy

The stabilization energy of the $N$-monomer cluster relative to infinitely separated monomers:

$$\boxed{
\Delta E_{\mathrm{bind}} \;=\;
\left(\, E_{\mathrm{M}} \;-\; \sum_{i=1}^{N} E_{\mathrm{mono\text{-}i}} \,\right)
\;\times\; 627.51 \;\;\text{kcal mol}^{-1}
}$$

**Average binding energy per hydrogen bond** ($n_{\mathrm{HB}}$ = number of detected HBs):

$$\boxed{
\overline{E}_{\mathrm{bind/HB}} \;=\; \frac{\Delta E_{\mathrm{bind}}}{n_{\mathrm{HB}}}
}$$

**Average MTA energy per hydrogen bond:**

$$\boxed{
\overline{E}_{\mathrm{MTA/HB}} \;=\; \frac{1}{n_{\mathrm{HB}}} \sum_{n=1}^{n_{\mathrm{HB}}} E_{\mathrm{MTA\text{-}HB}n}
}$$

---

### 5. MTA Error — Internal Consistency Check

The MTA error quantifies the deviation between the total binding energy and the sum of individual HBEs, serving as a self-consistency diagnostic:

$$\boxed{
\varepsilon_{\mathrm{MTA}} \;=\;
\left|\; |\Delta E_{\mathrm{bind}}| \;-\; \sum_{n=1}^{n_{\mathrm{HB}}} E_{\mathrm{MTA\text{-}HB}n} \;\right|
}$$

**Per-HB MTA error:**

$$\boxed{
\varepsilon_{\mathrm{MTA/HB}} \;=\; \frac{\varepsilon_{\mathrm{MTA}}}{n_{\mathrm{HB}}}
}$$

A small per-HB error (ideally $\varepsilon_{\mathrm{MTA/HB}} < 1\;\text{kcal mol}^{-1}$) confirms that the fragmentation scheme is self-consistent and that no hydrogen bonds have been missed or double-counted.

---

## Example

The repository includes a worked example: the cyclic hydrogen-bonded trimer **F1O2** (HF · 2H₂O), comprising 8 atoms and 3 monomers with 3 hydrogen bonds, computed at MP2/aug-cc-pVTZ with `scf=tight`.

**Input: `F1O2.xyz`**

```
8
F1O2
O    1.636242   0.060448   0.103696
H    1.108561  -0.755119   0.053012
H    2.377032  -0.081664  -0.498399
O   -0.852903   1.235815  -0.127148
H   -1.144618   1.756688   0.632186
H    0.122702   1.191890  -0.054566
F   -0.867229  -1.330843   0.007003
H   -1.042224  -0.390985  -0.036243
```

**Results (from `F1O2_summary.txt`):**

```
  Run Time          : 2026-03-11 16:36:24                 No. of monomers     : 3
  File              : F1O2.xyz                            No. of HB           : 3
  Mode              : calc                                Binding energy      : 19.891
  Level of Theory   : MP2/aug-cc-pVTZ  scf=tight          Avg per HB          : 6.630
  hb_min_dist       : 1.2                                 MTA_HBEs Sum        : 25.672
  hb_max_dist       : 2.3                                 Avg per MTA_HB      : 8.557
  hb_min_angle      : 120.0                               Total MTA error     : 5.781
  hb_max_angle      : 180.0                               Per-HB error        : 1.927
  neighbor_th       : 2.5                                 # All energies in kcal/mol

+-----+------------------+----------+----------+----------+-----------+---------------+---------------+---------------+
| HB# |     D-H...A      |  D-H(A)  | H...A(A) | A...D(A) |Angle(deg) |    MTA_HBE    |     Dimer     |     Coop      |
+-----+------------------+----------+----------+----------+-----------+---------------+---------------+---------------+
|  1  |    O1-H3...F7    |  0.973   |  2.058   |  2.866   |  139.144  |     6.309     |     3.418     |     2.890     |
|  2  |    O4-H5...O1    |  0.979   |  1.896   |  2.762   |  145.904  |     7.623     |     4.732     |     2.890     |
|  3  |    F7-H8...O4    |  0.957   |  1.640   |  2.570   |  162.827  |    11.740     |     8.850     |     2.890     |
+-----+------------------+----------+----------+----------+-----------+---------------+---------------+---------------+
```

All energies in kcal/mol. Distances in Ångströms.

HB3 (F–H···O, strongest at 11.7 kcal/mol) is the dominant interaction. All three HBs show equal cooperativity (2.890 kcal/mol), confirming the cyclic cooperative character of the trimer. The MTA error of 1.927 kcal/mol per HB is consistent with the basis set superposition inherent to the aug-cc-pVTZ fragmentation.

---

## Citations

If you use iHBEQuant in published work, please cite:

1. Patkar D; Ahirwar MB; Deshmukh MM *ChemPhysChem* **2022**, *23*, e202200476.
2. Patkar D; Ahirwar MB; Deshmukh MM *ChemPhysChem* **2022**, *23*, e202200143.
3. Patkar D; Ahirwar MB; Deshmukh MM *New J. Chem.* **2022**, *46*, 2368–2379.

---

## Author

**Deepak Patkar**  

For questions, bug reports, or feature suggestions, please open a [GitHub Issue](../../issues) or submit a pull request.

---

## License

This project is licensed under the [MIT License](LICENSE).
