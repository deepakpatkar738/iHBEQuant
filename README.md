<p align="center">
  <img src="iHBEQuant_logo.png" alt="iHBEQuant Logo" width="320"/>
</p>

<h1 align="center">Individual Hydrogen Bond Energy Quantifier</h1>

<p align="center">
  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.7%2B-blue.svg" alt="Python 3.7+"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green.svg" alt="License: MIT"></a>
  <img src="https://img.shields.io/badge/version-1.1-orange.svg" alt="Version">
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20Windows%20%7C%20macOS-lightgrey.svg" alt="Platform">
</p>

**iHBEQuant** is a Python tool that automates the quantification of individual hydrogen bond (HB) energies in molecular clusters using the **Molecular Tailoring Approach (MTA)**. Given a cluster geometry in XYZ format, it automatically detects hydrogen bonds, generates all required Gaussian input files, optionally submits and monitors calculations, and parses the results into a clean energy summary — all driven by a single plain-text configuration file.

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

**iHBEQuant** solves this by implementing the **Molecular Tailoring Approach (MTA)** for individual hydrogen bond energy (HBE) quantification, originally proposed by:

> Ahirwar MB; Gadre SR; Deshmukh MM *J. Phys. Chem. A* **2020**, *124*, 6699–6706. 

In MTA, the cluster is systematically fragmented — by removing the donor monomer, the acceptor monomer, or both — and the individual HBE is recovered from the total and fragment energies without any empirical parameters. iHBEQuant fully automates this fragmentation, Gaussian job management, energy parsing, and reporting workflow.

---

## Features

- **Automatic HB detection** — geometry-driven D–H···A search with fully configurable distance and angle thresholds
- **Fragment generation** — seven Gaussian input (`.gjf`) files generated per hydrogen bond: `cluster`, `rmD`, `rmA`, `rmDA`, `dimerDA`, `dimerD`, `dimerA`
- **Three run modes** — generate input files only (`gen`), auto-submit to Gaussian 16 and parse (`run`), or parse existing log files without resubmission (`cal`)
- **Cooperative energy** — decomposes each individual HBE into an isolated-dimer component and a cooperativity contribution
- **Binding energy and MTA error** — computes the total cluster binding energy and cross-validates it against the sum of individual HBEs
- **Multi-file batch processing** — accepts a space-separated list of `.xyz` files in `INPUT.cfg`; processes them sequentially
- **Supports DFT, MP2, CCSD(T)** — compatible with any single-point method available in Gaussian 16
- **Zero dependencies** — uses Python standard library only; no `pip install` required

---

## Requirements

- **Python 3.7+**
- **[Gaussian 16](https://gaussian.com/)** with `g16` accessible on the system PATH — required only for `run` and `cal` modes; `gen` mode works without Gaussian

---

## Installation

```bash
git clone https://github.com/deepakpatkar738/iHBEQuant.git
cd iHBEQuant
```

No compilation or package installation is needed. Run the script directly:

```bash
python iHBEQuant.py
```

---

## Quick Start

1. Place your cluster geometry as a standard XYZ file (e.g., `F1N2.xyz`) in the working directory.
2. Edit `INPUT.cfg` — set the XYZ filename, method, basis set, charge, multiplicity, and run mode.
3. Run:

```bash
python iHBEQuant.py
```

Results are written to a subdirectory named after your XYZ file (e.g., `F1N2/`).

---

## Usage

```bash
# Use default INPUT.cfg in the current directory
python iHBEQuant.py

# Specify a custom configuration file
python iHBEQuant.py -c myrun.cfg

# Verbose mode — prints the name of every .gjf file as it is written
python iHBEQuant.py -c myrun.cfg -v
```

### Command-Line Arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-c`, `--config` | `INPUT.cfg` | Path to the configuration file |
| `-v`, `--verbose` | off | Print all generated `.gjf` filenames to screen |

---

## Configuration

All run parameters are controlled through a plain-text `INPUT.cfg` file. Inline `#` comments are supported and stripped automatically. Sections can appear with or without blank lines between them.

```ini
[SYSTEM]
xyz_file  = F1N2.xyz        # space-separated list of .xyz files; processed one by one
nproc     = 8               # CPU cores for Gaussian
mem       = 8GB             # memory for Gaussian

[THEORY]
method       = MP2          # DFT functional, MP2, CCSD(T), etc.
basis        = aug-cc-pVTZ  # basis set
keywords     = scf=tight    # additional Gaussian route keywords (optional)
charge       = 0            # total charge of the full cluster (0 = neutral)
multiplicity = 1            # spin multiplicity (1 = singlet, 2 = doublet, ...)

[THRESHOLDS]
hb_min_dist  = 1.2          # H···A minimum distance (Å)
hb_max_dist  = 2.8          # H···A maximum distance (Å)
hb_min_angle = 100.0        # D–H···A minimum angle (degrees)
hb_max_angle = 180.0        # D–H···A maximum angle (degrees)
neighbor_th  = 2.5          # covalent neighbor cutoff (Å)

[RUN]
mode = gen
# gen  →  write .gjf input files only         (no Gaussian run)
# run  →  write .gjf files + submit to g16    (full calculation)
# cal  →  parse existing .log files only      (no new submission)
```

> **Important:** `charge` and `multiplicity` apply to the full cluster. All fragment calculations (rmD, rmA, dimerDA, etc.) inherit the same values. If your system contains ions or open-shell species, adjust these fields accordingly. iHBEQuant does **not** auto-detect charge or multiplicity from the XYZ file.

> **Tip:** Extra Gaussian route keywords (e.g., `scf=tight`, `int=ultrafine`) can be placed in the `keywords` field or appended directly to the `basis` field — both are handled correctly.

---

## Output

For each input XYZ file (e.g., `F1N2.xyz`), iHBEQuant creates a subdirectory `F1N2/` containing:

### Gaussian Input Files

| File | Contents |
|------|----------|
| `F1N2_cluster.gjf` | Full N-monomer cluster (single-point) |
| `F1N2_M1.gjf`, `_M2.gjf`, … | Individual monomers (for total binding energy) |
| `F1N2_HB1_rmD.gjf` | Cluster with donor monomer removed |
| `F1N2_HB1_rmA.gjf` | Cluster with acceptor monomer removed |
| `F1N2_HB1_rmDA.gjf` | Cluster with both donor and acceptor removed |
| `F1N2_HB1_dimerDA.gjf` | Donor + acceptor dimer |
| `F1N2_HB1_dimerD.gjf` | Donor monomer only |
| `F1N2_HB1_dimerA.gjf` | Acceptor monomer only |

Files `HB1_*` are replicated as `HB2_*`, `HB3_*`, … for every detected hydrogen bond.

### Summary and Coordinates

| File | Contents |
|------|----------|
| `F1N2_summary.txt` | Full energy table, run header, and per-HB calculation details |
| `F1N2_coords.xyz` | All fragment geometries in multi-block XYZ format for visualization |

---

## Methodology

iHBEQuant implements the **Molecular Tailoring Approach** for individual hydrogen bond energies as formulated by Ahirwar, Gadre, and Deshmukh (2020). For a hydrogen bond between donor monomer D and acceptor monomer A in an N-monomer cluster:

```
MTA_HBE  =  | E_cluster − ( E_rmD + E_rmA − E_rmDA ) |  ×  627.51 kcal/mol

Dimer    =  | E_dimerD + E_dimerA − E_dimerDA |  ×  627.51 kcal/mol

Coop     =  | MTA_HBE − Dimer |
```

**Symbol definitions:**

| Symbol | Meaning |
|--------|---------|
| `E_cluster` | Total energy of the full N-monomer cluster |
| `E_rmD` | Energy of the cluster with the donor monomer removed |
| `E_rmA` | Energy of the cluster with the acceptor monomer removed |
| `E_rmDA` | Energy of the cluster with both donor and acceptor removed |
| `E_dimerDA` | Energy of the isolated donor–acceptor dimer |
| `E_dimerD` | Energy of the isolated donor monomer |
| `E_dimerA` | Energy of the isolated acceptor monomer |

**MTA_HBE** captures the full in-cluster contribution of a hydrogen bond, including effects from all other monomers. **Dimer** is the pairwise interaction energy of the donor–acceptor pair in isolation. **Coop** is the difference between the two, quantifying how much the surrounding cluster environment strengthens or weakens the hydrogen bond — i.e., the cooperativity (Coop > 0) or anti-cooperativity (Coop < 0).

The **MTA error** provides an internal consistency check:

```
MTA error  =  | Binding energy | − Σ(MTA_HBE)
```

A small per-HB MTA error (ideally < 1 kcal/mol) confirms that the fragmentation scheme is self-consistent.

---

## Example

The repository includes a worked example: the hydrogen-bonded trimer **F1N2** (HF · 2NH₃), comprising 10 atoms and 3 monomers with 3 hydrogen bonds, computed at MP2/aug-cc-pVTZ.

**Input: `F1N2.xyz`**

```
10
N2H7F
F   -0.897087   1.409019   0.000010
H   -1.064941   0.446720   0.000039
N    1.854670  -0.137672  -0.000046
H    2.462588  -0.055488   0.815054
H    1.256831   0.693215  -0.000709
H    2.464218  -0.055691  -0.813946
N   -1.029638  -1.167465  -0.000014
H   -1.405133  -1.642027  -0.819856
H   -1.405101  -1.642103   0.819800
H   -0.009905  -1.289832  -0.000044
```

**Results (from `F1N2_summary.txt`):**

```
  File            : F1N2.xyz                   Mode : calc    Run Time : 2026-03-09 20:17:40
  Level of Theory : MP2/aug-cc-pVTZ scf=tight  No. of monomers: 3     No. of HB: 3
  Binding energy  : 22.216 kcal/mol            Avg per HB     : 7.405 kcal/mol
  MTA_HBEs Sum    : 26.598 kcal/mol
  Total MTA error : 4.382 kcal/mol             Per-HB error   : 1.461 kcal/mol

+-----+------------------+----------+----------+-----------+---------------+---------------+---------------+
| HB# |     D-H...A      |  D-H(A)  | H...A(A) |Angle(deg) |    MTA_HBE    |     Dimer     |     Coop      |
+-----+------------------+----------+----------+-----------+---------------+---------------+---------------+
|  1  |    F1-H2...N7    |  0.977   |  1.615   |  168.853  |    16.292     |    14.101     |     2.191     |
|  2  |    N3-H5...F1    |  1.024   |  2.270   |  144.119  |     5.048     |     2.857     |     2.191     |
|  3  |   N7-H10...N3    |  1.027   |  2.192   |  141.444  |     5.258     |     3.067     |     2.191     |
+-----+------------------+----------+----------+-----------+---------------+---------------+---------------+
```

All energies in kcal/mol. Distances in Ångströms.

HB1 (F–H···N, strongest at 16.3 kcal/mol) is the primary hydrogen bond; HB2 and HB3 are N–H···F and N–H···N contacts with significant cooperative enhancement (~2.2 kcal/mol each).

---

## Citations

If you use iHBEQuant in published work, please cite:

1. Ahirwar MB; Gadre SR; Deshmukh MM *J. Phys. Chem. A* **2020**, *124*, 6699–6706. 
2. Patkar D; Ahirwar MB; Deshmukh MM *ChemPhysChem* **2022**, *23*, e202200476.
3. Patkar D; Ahirwar MB; Deshmukh MM *ChemPhysChem* **2022**, *23*, e202200143.
4. Patkar D; Ahirwar MB; Deshmukh MM *New J. Chem.* **2022**, *46*, 2368–2379.
5. Patkar D; Ahirwar MB; Gadre SR; Deshmukh MM *J. Phys. Chem. A* **2021**, *125*, 8836–8845.

---

## Author

**Dr. Deepak Patkar**  

For questions, bug reports, or feature suggestions, please open a [GitHub Issue](../../issues) or submit a pull request.

---

## License

This project is licensed under the [MIT License](LICENSE).
