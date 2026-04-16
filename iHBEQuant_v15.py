#!/usr/bin/env python3
"""
iHBEQuant_F.py - Individual Hydrogen Bond Energy Quantifier
Version: 1.5 (fully configurable thresholds)
Author : Deepak Patkar

All thresholds are read from INPUT.cfg:
  covtol_cutoff, HB_acceptor, HB_donor, neigh_H_cutoff,
  hb_distance_min, hb_distance_max, hb_angle_min, hb_angle_max
"""

import os, sys, math, shutil, configparser, argparse, re, subprocess
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

G16    = "g16"
H2KCAL = 627.5094740631

LOGO = r"""
+---------------------------------------------------------------------------------------+
|                                                                                       |
|        _   ____  ____  ______   ________    ___                              _        |
|       (_) |_   ||   _||_   _ \ |_   __  | .'   `.                           / |_      |
|       __    | |__| |    | |_) |  | |_ \_|/  .-.  \  __   _   ,--.   _ .--. `| |-'     |
|      [  |   |  __  |    |  __'.  |  _| _ | |   | | [  | | | `'_\ : [ `.-. | | |       |
|       | |  _| |  | |_  _| |__) |_| |__/ |\  `-'  \_ | \_/ |,// | |, | | | | | |,      |
|      [___]|____||____||_______/|________| `.___.\__|'.__.'_/\'-;__/[___||__]\__/      |
|                                                                                       |
|      iHBEQuant : Individual Hydrogen Bond Energy Quantifier | version 1.5             |
|      Author    : Deepak Patkar                              | last update 2026        |
+---------------------------------------------------------------------------------------+
|  Cite us:                                                                             |
|  1. Patkar D; Ahirwar MB; Deshmukh MM ChemPhysChem 2022, 23, e202200476.              |
|  2. Patkar D; Ahirwar MB; Deshmukh MM ChemPhysChem 2022, 23, e202200143.              |
|  3. Patkar D; Ahirwar MB; Deshmukh MM New J. Chem. 2022, 46, 2368-2379.               |
+---------------------------------------------------------------------------------------+
"""

# ----------------------------------------------------------------------
# Configuration parser (reads all new thresholds)
# ----------------------------------------------------------------------
def load_cfg(path):
    if not os.path.isfile(path):
        sys.exit(f"ERROR: config '{path}' not found.")

    with open(path, encoding="latin-1") as f:
        raw_lines = f.readlines()

    # Remove comments and strip
    clean_lines = []
    for line in raw_lines:
        line = re.sub(r'\s*#.*', '', line).rstrip()
        if line:
            clean_lines.append(line)

    # --- Extract xyz_file list manually (supports brackets and multi-line) ---
    xyz_list = []
    in_system = False
    in_xyz = False
    bracket_depth = 0
    xyz_buffer = []

    for line in clean_lines:
        if line.startswith('[') and line.endswith(']'):
            in_system = (line.lower() == '[system]')
            in_xyz = False
            bracket_depth = 0
            continue
        if not in_system:
            continue
        if not in_xyz and line.startswith('xyz_file'):
            parts = line.split('=', 1)
            if len(parts) == 2:
                val = parts[1].strip()
                if val.startswith('['):
                    bracket_depth += 1
                    in_xyz = True
                    xyz_buffer = [val]
                else:
                    xyz_list = re.findall(r'[\w.-]+\.xyz', val)
                    break
        elif in_xyz:
            xyz_buffer.append(line)
            bracket_depth += line.count('[') - line.count(']')
            if bracket_depth == 0:
                full_value = ' '.join(xyz_buffer)
                xyz_list = re.findall(r'[\w.-]+\.xyz', full_value)
                break

    # Fallback: single-line without brackets
    if not xyz_list:
        for line in clean_lines:
            if line.startswith('xyz_file'):
                parts = line.split('=', 1)
                if len(parts) == 2:
                    xyz_list = re.findall(r'[\w.-]+\.xyz', parts[1])
                break

    # --- Prepare a clean config string for configparser (replace xyz_file with dummy) ---
    config_lines = []
    for line in clean_lines:
        if line.startswith('xyz_file'):
            config_lines.append('xyz_file = dummy.xyz')
        else:
            config_lines.append(line)
    config_str = "\n".join(config_lines)
    c = configparser.ConfigParser()
    c.read_string(config_str)

    MODE = {"gen":"generate","generate":"generate","run":"run","submit":"run",
            "cal":"calc","calc":"calc"}
    chk_opt = c.get("SYSTEM", "chk", fallback="N").strip().upper()
    frq_opt = c.get("THEORY", "frq", fallback="N").strip().upper()

    # Read thresholds with new names, fallback to old names for compatibility
    # covtol_cutoff
    covtol = c.getfloat("THRESHOLDS", "covtol_cutoff", fallback=None)
    if covtol is None:
        covtol = c.getfloat("THRESHOLDS", "cov_tol", fallback=0.5)   # old name fallback

    # HB_acceptor list
    acc_str = c.get("THRESHOLDS", "HB_acceptor", fallback="F,O,N,S")
    HB_acceptor = set(re.split(r'[,\s]+', acc_str.strip().upper()))

    # HB_donor list
    don_str = c.get("THRESHOLDS", "HB_donor", fallback="F,O,N,S,C,Cl")
    HB_donor = set(re.split(r'[,\s]+', don_str.strip().upper()))

    neigh_H = c.getfloat("THRESHOLDS", "neigh_H_cutoff", fallback=2.5)

    # Distance min/max
    d_min = c.getfloat("THRESHOLDS", "hb_distance_min", fallback=None)
    if d_min is None:
        d_min = c.getfloat("THRESHOLDS", "hb_min_dist", fallback=1.2)
    d_max = c.getfloat("THRESHOLDS", "hb_distance_max", fallback=None)
    if d_max is None:
        d_max = c.getfloat("THRESHOLDS", "hb_max_dist", fallback=2.5)

    # Angle min/max
    a_min = c.getfloat("THRESHOLDS", "hb_angle_min", fallback=None)
    if a_min is None:
        a_min = c.getfloat("THRESHOLDS", "hb_min_angle", fallback=100.0)
    a_max = c.getfloat("THRESHOLDS", "hb_angle_max", fallback=None)
    if a_max is None:
        a_max = c.getfloat("THRESHOLDS", "hb_max_angle", fallback=180.0)

    return {
        "xyz_list" : xyz_list,
        "nproc"    : c.getint("SYSTEM", "nproc", fallback=8),
        "mem"      : c.get("SYSTEM", "mem", fallback="4GB").strip(),
        "method"   : c.get("THEORY", "method", fallback="B3LYP").strip(),
        "basis"    : c.get("THEORY", "basis", fallback="6-311+G(d,p)").strip(),
        "keywords" : c.get("THEORY", "keywords", fallback="").strip(),
        "charge"   : c.getint("THEORY", "charge", fallback=0),
        "mult"     : c.getint("THEORY", "multiplicity", fallback=1),
        "covtol_cutoff" : covtol,
        "HB_acceptor"   : HB_acceptor,
        "HB_donor"      : HB_donor,
        "neigh_H_cutoff": neigh_H,
        "d_min"    : d_min,
        "d_max"    : d_max,
        "a_min"    : a_min,
        "a_max"    : a_max,
        "mode"     : MODE.get(c.get("RUN", "mode", fallback="gen").strip().lower(), "generate"),
        "njob"     : c.getint("RUN", "njob", fallback=1),
        "chk"      : chk_opt,
        "frq"      : frq_opt,
    }

# ----------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------
def load_xyz(path):
    raw = open(path).readlines()
    def _parse(lines):
        out = []
        for l in lines:
            p = l.split()
            if len(p) < 4 or p[0].lstrip("+-").replace(".", "").isdigit():
                continue
            try:
                x, y, z = float(p[1]), float(p[2]), float(p[3])
            except ValueError:
                continue
            out.append({"s": p[0], "x": x, "y": y, "z": z,
                        "gjf": f"{p[0]:<4}{x:>14.8f}{y:>14.8f}{z:>14.8f}\n"})
        return out
    try:
        n = int(raw[0].strip())
        return _parse(raw[2:2+n])
    except (ValueError, IndexError):
        return _parse(raw)

def load_gjf_atoms(gjf_path):
    """Read Gaussian .gjf and return list of atom dicts (same format as load_xyz)."""
    if not os.path.isfile(gjf_path):
        sys.exit(f"ERROR: cannot read parent gjf '{gjf_path}' - file missing.")
    atoms = []
    blank_count = 0
    charge_mult_seen = False
    with open(gjf_path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('%') or stripped.startswith('#'):
                continue
            if stripped == '':
                blank_count += 1
                continue
            if blank_count >= 2 and not charge_mult_seen:
                parts = stripped.split()
                if len(parts) == 2 and all(p.lstrip('-').isdigit() for p in parts):
                    charge_mult_seen = True
                    continue
            if charge_mult_seen:
                parts = stripped.split()
                if len(parts) < 4:
                    break
                try:
                    x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                except ValueError:
                    break
                sym = parts[0]
                atoms.append({"s": sym, "x": x, "y": y, "z": z,
                              "gjf": f"{sym:<4}{x:>14.8f}{y:>14.8f}{z:>14.8f}\n"})
    if not atoms:
        sys.exit(f"ERROR: no atoms parsed from '{gjf_path}'.")
    return atoms

def dist(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

def lbl(i, atoms):
    return f"{atoms[i]['s']}{i+1}"

def mono_label(mi, monomers, atoms):
    return "-".join(lbl(a, atoms) for a in monomers[mi])

def bonded(a, b, cfg):
    """Use covalent radii sum + covtol_cutoff from config."""
    r_cov = 0.0
    for sym in (a["s"].upper(), b["s"].upper()):
        # simple covalent radii lookup (common elements)
        rad = {"H":0.23, "C":0.77, "N":0.75, "O":0.70, "F":0.71,
               "S":1.05, "P":1.06, "CL":0.99, "SI":1.11}.get(sym, 0.90)
        r_cov += rad
    tol = cfg["covtol_cutoff"]
    return abs(dist(a, b) - r_cov) <= tol

def hb_angle(donor, h, acc):
    def v(o, t):
        return [t["x"]-o["x"], t["y"]-o["y"], t["z"]-o["z"]]
    v1, v2 = v(h, donor), v(h, acc)
    dot = sum(i*j for i,j in zip(v1, v2))
    nm = lambda u: math.sqrt(sum(i*i for i in u))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot/(nm(v1)*nm(v2)+1e-12)))))

def find_monomers(atoms, cfg):
    visited, monomers = set(), []
    def fill(start):
        stack, m = [start], []
        while stack:
            i = stack.pop()
            if i in visited:
                continue
            visited.add(i)
            m.append(i)
            for j in range(len(atoms)):
                if j not in visited and bonded(atoms[i], atoms[j], cfg):
                    stack.append(j)
        return sorted(m)
    for i in range(len(atoms)):
        if i not in visited:
            monomers.append(fill(i))
    return monomers

def find_hbonds(atoms, monomers, cfg):
    mono_of = {a: mi for mi, m in enumerate(monomers) for a in m}
    hbonds = []
    # Precompute acceptor set and donor set from config
    acceptor_set = cfg["HB_acceptor"]
    donor_set    = cfg["HB_donor"]

    for hi in range(len(atoms)):
        if atoms[hi]["s"].upper() != "H":
            continue
        # Find donor (heavy atom covalently bonded to this H)
        donor = None
        for d in range(len(atoms)):
            if d != hi and atoms[d]["s"].upper() != "H" and bonded(atoms[hi], atoms[d], cfg):
                # Check if donor element is allowed
                if atoms[d]["s"].upper() in donor_set:
                    donor = d
                    break
        if donor is None:
            continue

        # For each possible acceptor
        for ai in range(len(atoms)):
            if ai in (hi, donor):
                continue
            elem = atoms[ai]["s"].upper()
            if elem not in acceptor_set:
                continue
            if mono_of[ai] == mono_of[hi]:
                continue

            # Quick distance filter: H...A must be <= neigh_H_cutoff (first sphere)
            ha = dist(atoms[hi], atoms[ai])
            if ha > cfg["neigh_H_cutoff"]:
                continue

            # Now apply the strict HB distance and angle criteria
            if not (cfg["d_min"] < ha < cfg["d_max"]):
                continue
            ang = hb_angle(atoms[donor], atoms[hi], atoms[ai])
            if not (cfg["a_min"] <= ang <= cfg["a_max"]):
                continue

            hbonds.append({
                "di": donor, "hi": hi, "ai": ai,
                "dm": mono_of[donor], "am": mono_of[ai],
                "dh_dist": round(dist(atoms[donor], atoms[hi]), 3),
                "ha_dist": round(ha, 3),
                "angle": round(ang, 3),
                "da_dist": round(dist(atoms[donor], atoms[ai]), 3)
            })
    return hbonds

# ----------------------------------------------------------------------
# Gaussian interface (unchanged except route_line)
# ----------------------------------------------------------------------
def route_line(cfg, freq=False):
    kw = f" {cfg['keywords'].strip()}" if cfg.get("keywords", "").strip() else ""
    freq_str = " freq" if freq else ""
    return f"#P {cfg['method']}/{cfg['basis']} SP nosymm{freq_str}{kw}"

def write_gjf(path, idx, atoms, cfg, title, freq=False):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"%chk={os.path.splitext(path)[0]}.chk\n")
        f.write(f"%nproc={cfg['nproc']}\n%mem={cfg['mem']}\n")
        f.write(f"{route_line(cfg, freq)}\n\n{title}\n\n{cfg['charge']} {cfg['mult']}\n")
        for i in idx:
            f.write(atoms[i]["gjf"])
        f.write("\n\n")

def write_gjf_from_gjf(path, parent_atoms, local_indices, cfg, title, freq=False):
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"%chk={os.path.splitext(path)[0]}.chk\n")
        f.write(f"%nproc={cfg['nproc']}\n%mem={cfg['mem']}\n")
        f.write(f"{route_line(cfg, freq)}\n\n{title}\n\n{cfg['charge']} {cfg['mult']}\n")
        for i in local_indices:
            f.write(parent_atoms[i]["gjf"])
        f.write("\n\n")

def log_path(gjf):
    return gjf.replace(".gjf", ".log")

def normal_end(log):
    return os.path.isfile(log) and any("Normal termination" in l for l in open(log))

def run_g16(gjf):
    if normal_end(log_path(gjf)):
        return True
    print(f"  Running g16: {os.path.basename(gjf)}")
    ret = subprocess.run([G16, gjf], capture_output=False)
    return ret.returncode == 0

def run_g16_parallel(gjf_list, max_jobs):
    if not gjf_list:
        return
    print(f"\n  Running {len(gjf_list)} Gaussian jobs with concurrency = {max_jobs} ...")
    with ThreadPoolExecutor(max_workers=max_jobs) as executor:
        futures = {executor.submit(run_g16, gjf): gjf for gjf in gjf_list}
        for future in as_completed(futures):
            gjf = futures[future]
            try:
                if not future.result():
                    print(f"  WARNING: g16 non-zero exit for {gjf}")
            except Exception as e:
                print(f"  ERROR running {gjf}: {e}")

def parse_energy(log):
    e = None
    for l in open(log):
        if "SCF Done" in l:
            e = float(l.split()[4])
        if "EUMP2" in l:
            e = float(l.split()[5].replace("D", "E"))
        if "E(Corr)=" in l:
            e = float(l.split()[4])
        if "CCSD(T)= " in l:
            e = float(l.split()[1].replace("D", "E"))
    return e

def parse_dipole(log):
    if not os.path.isfile(log):
        return None
    dipole = None
    try:
        with open(log) as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if "Dipole moment (field-independent basis, Debye)" in line:
                if i + 1 < len(lines):
                    val_line = lines[i + 1]
                    m = re.search(r'Tot=\s*([-\d.]+)', val_line)
                    if m:
                        dipole = float(m.group(1))
    except Exception:
        pass
    return dipole

def get_energy(gjf, mode, run_later_list=None):
    if mode == "generate":
        return None
    if mode == "run":
        if normal_end(log_path(gjf)):
            return parse_energy(log_path(gjf))
        else:
            if run_later_list is not None:
                run_later_list.append(gjf)
            return None
    # mode == "calc"
    lp = log_path(gjf)
    if normal_end(lp):
        return parse_energy(lp)
    print(f"  WARNING: abnormal/missing log - {gjf}")
    return None

def make_fchk(gjf):
    chk_file = gjf.replace(".gjf", ".chk")
    fchk_file = gjf.replace(".gjf", ".fchk")
    if os.path.isfile(chk_file) and not os.path.isfile(fchk_file):
        print(f"  Creating .fchk: {os.path.basename(fchk_file)}")
        subprocess.run(["formchk", chk_file, fchk_file], capture_output=True)

# ----------------------------------------------------------------------
# File naming (unchanged)
# ----------------------------------------------------------------------
def parent_gjf(out_dir, prefix):
    return f"{out_dir}/{prefix}_M.gjf"

def mono_gjf(out_dir, prefix, mi):
    return f"{out_dir}/{prefix}_mono_M{mi+1}.gjf"

def frag_gjf(out_dir, prefix, mi):
    return f"{out_dir}/{prefix}_frag_M{mi+1}.gjf"

def overlap_gjf(out_dir, prefix, mi, mj):
    i, j = sorted([mi+1, mj+1])
    return f"{out_dir}/{prefix}_frag_M{i}M{j}.gjf"

def dimer_gjf(out_dir, prefix, mi, mj):
    i, j = sorted([mi+1, mj+1])
    return f"{out_dir}/{prefix}_dimer_M{i}M{j}.gjf"

# ----------------------------------------------------------------------
# Fragment generation (all read from _M.gjf) - unchanged
# ----------------------------------------------------------------------
def generate_primary_fragments(prefix, out_dir, parent_atoms, monomers, cfg, atoms):
    all_idx = [i for m in monomers for i in m]
    pos_in_parent = {xyz_i: k for k, xyz_i in enumerate(all_idx)}
    frag_files = {}
    for mi, mon_idx in enumerate(monomers):
        frag_xyz_idx = [i for mj, m in enumerate(monomers) for i in m if mj != mi]
        frag_parent_idx = [pos_in_parent[i] for i in frag_xyz_idx]
        gjf = frag_gjf(out_dir, prefix, mi)
        if not os.path.isfile(gjf):
            title = f"{prefix} primary frag | removed monomer M{mi+1} ({mono_label(mi, monomers, atoms)})"
            write_gjf_from_gjf(gjf, parent_atoms, frag_parent_idx, cfg, title, freq=False)
        frag_files[mi] = gjf
    return frag_files

def generate_pair_fragments(prefix, out_dir, parent_atoms, monomers, cfg, pair_set, atoms):
    all_idx = [i for m in monomers for i in m]
    pos_in_parent = {xyz_i: k for k, xyz_i in enumerate(all_idx)}
    overlap_files = {}
    dimer_files   = {}
    for mi, mj in pair_set:
        ov_xyz_idx = [i for mk, m in enumerate(monomers) for i in m if mk not in (mi, mj)]
        ov_parent_idx = [pos_in_parent[i] for i in ov_xyz_idx]
        ov_gjf = overlap_gjf(out_dir, prefix, mi, mj)
        if not os.path.isfile(ov_gjf):
            title = f"{prefix} overlap frag | removed M{mi+1} and M{mj+1}"
            write_gjf_from_gjf(ov_gjf, parent_atoms, ov_parent_idx, cfg, title, freq=False)
        overlap_files[(mi, mj)] = ov_gjf

        dim_xyz_idx = monomers[mi] + monomers[mj]
        dim_parent_idx = [pos_in_parent[i] for i in dim_xyz_idx]
        dim_gjf = dimer_gjf(out_dir, prefix, mi, mj)
        if not os.path.isfile(dim_gjf):
            title = f"{prefix} dimer M{mi+1}M{mj+1}"
            write_gjf_from_gjf(dim_gjf, parent_atoms, dim_parent_idx, cfg, title, freq=False)
        dimer_files[(mi, mj)] = dim_gjf
    return overlap_files, dimer_files

def generate_monomers(prefix, out_dir, parent_atoms, monomers, cfg, atoms):
    all_idx = [i for m in monomers for i in m]
    pos_in_parent = {xyz_i: k for k, xyz_i in enumerate(all_idx)}
    mono_gjfs = []
    for mi, m in enumerate(monomers):
        gjf = mono_gjf(out_dir, prefix, mi)
        mono_parent_idx = [pos_in_parent[i] for i in m]
        write_gjf_from_gjf(gjf, parent_atoms, mono_parent_idx, cfg,
                            f"{prefix} monomer M{mi+1} | {mono_label(mi, monomers, atoms)}",
                            freq=False)
        mono_gjfs.append(gjf)
    return mono_gjfs

# ----------------------------------------------------------------------
# Energy calculation for a single HB (unchanged)
# ----------------------------------------------------------------------
def compute_hb_energy(hb, parent_e, mono_es, frag_es, overlap_es, dimer_es):
    mi = hb["dm"]
    mj = hb["am"]
    if None not in (frag_es[mi], frag_es[mj], overlap_es[(mi, mj)], parent_e):
        mta = (frag_es[mi] + frag_es[mj] - overlap_es[(mi, mj)]) - parent_e
        mta_kcal = round(mta * H2KCAL, 3)
    else:
        mta_kcal = None
    if None not in (mono_es[mi], mono_es[mj], dimer_es[(mi, mj)]):
        dimer = (mono_es[mi] + mono_es[mj] - dimer_es[(mi, mj)])
        dimer_kcal = round(dimer * H2KCAL, 3)
    else:
        dimer_kcal = None
    coop_kcal = round(mta_kcal - dimer_kcal, 3) if (mta_kcal is not None and dimer_kcal is not None) else None
    return mta_kcal, dimer_kcal, coop_kcal

# ----------------------------------------------------------------------
# Output formatting (modified: indices from .gjf order)
# ----------------------------------------------------------------------
def fmt(v, dec=3):
    return f"{v:.{dec}f}" if v is not None else "pending"

def xyz_block_from_gjf(gjf_path):
    if not os.path.isfile(gjf_path):
        return f"  [file not found: {gjf_path}]\n"
    gjf_atoms = load_gjf_atoms(gjf_path)
    title_from_gjf = os.path.basename(gjf_path)
    blank_count = 0
    with open(gjf_path, encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('%') or stripped.startswith('#'):
                continue
            if stripped == '':
                blank_count += 1
                continue
            if blank_count == 1:
                title_from_gjf = stripped
                break
    lines = [str(len(gjf_atoms)), title_from_gjf]
    for a in gjf_atoms:
        lines.append(f"{a['s']:<4}{a['x']:>14.8f}{a['y']:>14.8f}{a['z']:>14.8f}")
    lines.append("")
    return "\n".join(lines)

def dipole_tag(gjf_path):
    lp = log_path(gjf_path)
    d = parse_dipole(lp)
    if d is not None:
        return f"[Dipole]  {d:.4f} Debye"
    return "[Dipole]  pending"

def write_summary(path, records, atoms, monomers, mono_e, cluster_e, cfg, xyz, mode,
                  run_time, prefix, out_dir):
    n_hb   = len(records)
    n_mono = len(monomers)

    # --- Build mapping from original atom index to position in parent .gjf file ---
    all_idx = [i for m in monomers for i in m]   # order in which atoms appear in _M.gjf
    gjf_order = {orig: pos for pos, orig in enumerate(all_idx)}  # orig -> 0-based position

    def get_gjf_label(orig_idx):
        """Return label like 'F1' where number is 1-based position in parent .gjf file."""
        gjf_pos = gjf_order[orig_idx] + 1
        return f"{atoms[orig_idx]['s']}{gjf_pos}"

    if None not in [cluster_e] + mono_e and mode in ("run", "calc"):
        be_kcal = (cluster_e - sum(mono_e)) * H2KCAL
        avg_be  = be_kcal / n_hb if n_hb else None
    else:
        be_kcal = avg_be = None
    hbe_vals = [r["res"]["hbe"] for r in records if r["res"]["hbe"] is not None]
    sum_hbe  = round(sum(hbe_vals), 3) if hbe_vals else None
    if None not in (be_kcal, sum_hbe):
        total_err  = round(abs(abs(be_kcal) - sum_hbe), 3)
        per_hb_err = round(total_err / n_hb, 3) if n_hb else None
    else:
        total_err = per_hb_err = None

    lvl    = f"{cfg['method']}/{cfg['basis']}"
    kw_str = cfg.get("keywords", "").strip()
    if kw_str:
        lvl = f"{lvl}  {kw_str}"
    avg_mta = round(sum_hbe / n_hb, 3) if (sum_hbe is not None and n_hb) else None

    def left_right(left_label, left_val, right_label, right_val):
        left  = f"  {left_label:<18}: {str(left_val)}"
        right = f"  {right_label:<20}: {str(right_val)}" if right_label else ""
        return left.ljust(52) + right

    lines = [LOGO,
        left_right("Run Time",       run_time, "No. of monomers",         n_mono),
        left_right("File",           xyz,      "No. of HB",               n_hb),
        left_right("Mode",           mode,     "Binding energy",           fmt(abs(be_kcal) if be_kcal else None)),
        left_right("Level of Theory",lvl,      "Avg per HB",               fmt(abs(avg_be) if avg_be else None)),
        left_right("covtol_cutoff",  cfg["covtol_cutoff"], "MTA_HBEs Sum",         fmt(sum_hbe)),
        left_right("HB_acceptor",    ", ".join(sorted(cfg["HB_acceptor"])), "Avg per MTA_HB",       fmt(avg_mta)),
        left_right("HB_donor",       ", ".join(sorted(cfg["HB_donor"])), "Total MTA error",      fmt(total_err)),
        left_right("neigh_H_cutoff", cfg["neigh_H_cutoff"], "Per-HB error",         fmt(per_hb_err)),
        left_right("hb_distance",    f"{cfg['d_min']} - {cfg['d_max']}", "# All energies in kcal/mol", ""),
        left_right("hb_angle",       f"{cfg['a_min']} - {cfg['a_max']}", "", ""),
        ""]

    table_rows = []
    for r in records:
        hb, res = r["hb"], r["res"]
        # Use get_gjf_label for atom indices in HB type
        hb_type     = f"HB{r['n']} [{get_gjf_label(hb['di'])}-{get_gjf_label(hb['hi'])}...{get_gjf_label(hb['ai'])}]"
        donor_atoms = [get_gjf_label(i) for i in monomers[hb['dm']]]
        donor       = f"M{hb['dm']+1} [{' '.join(donor_atoms)}]"
        acc_atoms   = [get_gjf_label(i) for i in monomers[hb['am']]]
        acceptor    = f"M{hb['am']+1} [{' '.join(acc_atoms)}]"
        row = {
            "HB Type":    hb_type,
            "Donor":      donor,
            "Acceptor":   acceptor,
            "r(D-H)":     fmt(hb['dh_dist'], 3),
            "r(H...A)":   fmt(hb['ha_dist'], 3),
            "r(D...A)":   fmt(hb['da_dist'], 3),
            "a(D-H...A)": fmt(hb['angle'], 2),
            "HB_MTA":     fmt(res['hbe'],   3),
            "HB_Dimer":   fmt(res['dimer'], 3),
            "HB_Coop":    fmt(res['coop'],  3),
        }
        table_rows.append(row)

    if table_rows:
        col_names = list(table_rows[0].keys())
        col_width = {}
        for col in col_names:
            max_len = len(col)
            for row in table_rows:
                val = str(row[col])
                if len(val) > max_len:
                    max_len = len(val)
            col_width[col] = max_len + 2
        header = "  ".join(f"{col:<{col_width[col]}}" for col in col_names)
        lines.append(header)
        for row in table_rows:
            line = "  ".join(f"{str(row[col]):<{col_width[col]}}" for col in col_names)
            lines.append(line)
        lines.append("")

    lines += ["  CALCULATION DETAILS", "",
              f"  {xyz:<24} No of Atoms: {len(atoms):<8} No of Monomers: {n_mono:<8} HB detected: {n_hb}", ""]

    for r in records:
        hb, res = r["hb"], r["res"]
        mi, mj  = hb["dm"], hb["am"]
        tag     = f"HB{r['n']}"
        # Use gjf label for detailed per-HB line as well
        label   = f"{get_gjf_label(hb['di'])}-{get_gjf_label(hb['hi'])}...{get_gjf_label(hb['ai'])}"
        T       = "-"*95
        lines  += [f"  {tag}: {label}", T]

        def get_e(gjf):
            if mode == "generate":
                return None
            return get_energy(gjf, "calc")

        p_gjf_path  = parent_gjf(out_dir, prefix)
        fmi_gjf     = frag_gjf(out_dir, prefix, mi)
        fmj_gjf     = frag_gjf(out_dir, prefix, mj)
        ov_gjf_path = overlap_gjf(out_dir, prefix, mi, mj)
        dim_gjf_path= dimer_gjf(out_dir, prefix, mi, mj)
        mmi_gjf     = mono_gjf(out_dir, prefix, mi)
        mmj_gjf     = mono_gjf(out_dir, prefix, mj)

        parent_e_loc  = get_e(p_gjf_path)
        frag_mi_e     = get_e(fmi_gjf)
        frag_mj_e     = get_e(fmj_gjf)
        overlap_e_loc = get_e(ov_gjf_path)
        dimer_e_loc   = get_e(dim_gjf_path)
        mono_mi_e     = get_e(mmi_gjf)
        mono_mj_e     = get_e(mmj_gjf)

        def erow(ll, lv, rl=None, rv=None):
            left  = f"    {ll:<12}  {fmt(lv,6)}"
            right = f"    {rl:<12}  {fmt(rv,6)}" if rl else "                              # Hartree"
            lines.append(f"{left:<42}{right}")

        erow("parent_M",          parent_e_loc,  "dimer_M1M2",    dimer_e_loc)
        erow(f"frag_M{mi+1}",     frag_mi_e,     f"mono_M{mi+1}", mono_mi_e)
        erow(f"frag_M{mj+1}",     frag_mj_e,     f"mono_M{mj+1}", mono_mj_e)
        erow(f"frag_M{mi+1}M{mj+1}", overlap_e_loc)
        lines += [T,
                  f"    MTA_{tag}  : {fmt(res['hbe'])} kcal/mol   |   Dimer : {fmt(res['dimer'])} kcal/mol   |   Coop  : {fmt(res['coop'])} kcal/mol",
                  ""]

    lines += ["", "  COORDINATES", "  " + "="*91, ""]

    p_gjf_path  = parent_gjf(out_dir, prefix)
    all_count   = sum(len(m) for m in monomers)
    dp_parent   = dipole_tag(p_gjf_path)
    lines.append(f"  # Parent molecule : {prefix}_M  ({all_count} atoms, {n_mono} monomers)  {dp_parent}")
    lines.append(xyz_block_from_gjf(p_gjf_path))

    for mi in range(n_mono):
        mgi        = mono_gjf(out_dir, prefix, mi)
        dp_mono    = dipole_tag(mgi)
        n_atoms_mi = len(monomers[mi])
        lines.append(f"  # Monomer M{mi+1} : {mono_label(mi, monomers, atoms)}  ({n_atoms_mi} atoms)  {dp_mono}")
        lines.append(xyz_block_from_gjf(mgi))

    for mi in range(n_mono):
        fgi      = frag_gjf(out_dir, prefix, mi)
        dp_frag  = dipole_tag(fgi)
        lines.append(f"  # Primary fragment (remove M{mi+1}) : {os.path.basename(fgi)}  {dp_frag}")
        lines.append(xyz_block_from_gjf(fgi))

    seen_pairs = set()
    for r in records:
        mi, mj = r["hb"]["dm"], r["hb"]["am"]
        if (mi, mj) in seen_pairs or (mj, mi) in seen_pairs:
            continue
        seen_pairs.add((mi, mj))

        ov_gjf_path  = overlap_gjf(out_dir, prefix, mi, mj)
        dp_ov        = dipole_tag(ov_gjf_path)
        lines.append(f"  # Overlap fragment (remove M{mi+1} and M{mj+1}) : {os.path.basename(ov_gjf_path)}  {dp_ov}")
        lines.append(xyz_block_from_gjf(ov_gjf_path))

        dim_gjf_path = dimer_gjf(out_dir, prefix, mi, mj)
        dp_dim       = dipole_tag(dim_gjf_path)
        lines.append(f"  # Dimer (keep M{mi+1} and M{mj+1}) : {os.path.basename(dim_gjf_path)}  {dp_dim}")
        lines.append(xyz_block_from_gjf(dim_gjf_path))

    text = "\n".join(lines)
    print(text)
    open(path, "w", encoding="utf-8").write(text)

# ----------------------------------------------------------------------
# Main driver
# ----------------------------------------------------------------------
def run_one(xyz, cfg, verbose):
    mode     = cfg["mode"]
    run_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not os.path.isfile(xyz):
        print(f"  ERROR: '{xyz}' not found - skipping.\n")
        return

    prefix  = os.path.splitext(os.path.basename(xyz))[0]
    out_dir = prefix
    os.makedirs(out_dir, exist_ok=True)

    parent_xyz = f"{out_dir}/{prefix}_M.xyz"
    shutil.copy2(xyz, parent_xyz)
    print(f"  Parent molecule  : {parent_xyz}")

    atoms    = load_xyz(xyz)
    monomers = find_monomers(atoms, cfg)
    hbonds   = find_hbonds(atoms, monomers, cfg)

    print(f"  Monomers detected: {len(monomers)}")
    print(f"  H-bonds detected : {len(hbonds)}")

    if not hbonds:
        print(f"  No HB found in {xyz}. Adjust thresholds in INPUT.cfg.\n")
        return

    p_gjf   = parent_gjf(out_dir, prefix)
    all_idx = [i for m in monomers for i in m]
    do_freq = (cfg["frq"] == "Y")
    write_gjf(p_gjf, all_idx, atoms, cfg,
              f"{prefix} parent molecule | {len(atoms)} atoms {len(monomers)} monomers",
              freq=do_freq)

    parent_atoms = load_gjf_atoms(p_gjf)

    mono_gjfs = generate_monomers(prefix, out_dir, parent_atoms, monomers, cfg, atoms)
    frag_files = generate_primary_fragments(prefix, out_dir, parent_atoms, monomers, cfg, atoms)

    pair_set = set()
    for hb in hbonds:
        mi, mj = hb["dm"], hb["am"]
        if mi != mj:
            pair_set.add((mi, mj))

    overlap_files, dimer_files = generate_pair_fragments(
        prefix, out_dir, parent_atoms, monomers, cfg, pair_set, atoms)

    run_later_list = [] if mode == "run" else None
    energy_cache   = {}

    def get_cached_energy(gjf):
        if gjf in energy_cache:
            return energy_cache[gjf]
        e = get_energy(gjf, mode, run_later_list)
        energy_cache[gjf] = e
        return e

    parent_e  = get_cached_energy(p_gjf)
    mono_e    = [get_cached_energy(gjf) for gjf in mono_gjfs]
    frag_e    = {mi: get_cached_energy(gjf) for mi, gjf in frag_files.items()}
    overlap_e = {pair: get_cached_energy(gjf) for pair, gjf in overlap_files.items()}
    dimer_e   = {pair: get_cached_energy(gjf) for pair, gjf in dimer_files.items()}

    records = []
    for n, hb in enumerate(hbonds, 1):
        mi, mj = hb["dm"], hb["am"]
        pair   = (mi, mj) if (mi, mj) in overlap_e else (mj, mi)
        mta_kcal, dimer_kcal, coop_kcal = compute_hb_energy(
            hb, parent_e, mono_e, frag_e, overlap_e, dimer_e)
        records.append({
            "n":   n,
            "hb":  hb,
            "res": {"hbe": mta_kcal, "dimer": dimer_kcal, "coop": coop_kcal},
            "raw": {}
        })
        if verbose:
            print(f"    HB{n}: MTA={mta_kcal}, Dimer={dimer_kcal}, Coop={coop_kcal}")

    if mode == "run" and run_later_list:
        gjf_set = set(run_later_list)
        run_g16_parallel(list(gjf_set), cfg["njob"])
        energy_cache.clear()
        parent_e  = get_cached_energy(p_gjf)
        mono_e    = [get_cached_energy(gjf) for gjf in mono_gjfs]
        frag_e    = {mi: get_cached_energy(gjf) for mi, gjf in frag_files.items()}
        overlap_e = {pair: get_cached_energy(gjf) for pair, gjf in overlap_files.items()}
        dimer_e   = {pair: get_cached_energy(gjf) for pair, gjf in dimer_files.items()}
        for r in records:
            hb = r["hb"]
            mi, mj = hb["dm"], hb["am"]
            mta, dimer, coop = compute_hb_energy(hb, parent_e, mono_e, frag_e, overlap_e, dimer_e)
            r["res"]["hbe"]   = mta
            r["res"]["dimer"] = dimer
            r["res"]["coop"]  = coop

    if mode in ("run", "calc"):
        if cfg["chk"] == "Y":
            all_gjfs = ([p_gjf] + mono_gjfs +
                        list(frag_files.values()) +
                        list(overlap_files.values()) +
                        list(dimer_files.values()))
            for gjf in all_gjfs:
                make_fchk(gjf)
        else:
            make_fchk(p_gjf)

    summary = f"{out_dir}/{prefix}_summary.txt"
    write_summary(summary, records, atoms, monomers, mono_e, parent_e,
                  cfg, xyz, mode, run_time, prefix, out_dir)
    print(f"\n  Summary -> {summary}  (coordinates included)\n")

def main():
    ap = argparse.ArgumentParser(description="iHBEQuant: Individual HB Energy Quantifier")
    ap.add_argument("-c", "--config", default="INPUT.cfg", help="config file [INPUT.cfg]")
    ap.add_argument("-v", "--verbose", action="store_true", help="list gjf files written")
    args = ap.parse_args()
    cfg  = load_cfg(args.config)
    if cfg["mode"] not in ("generate", "run", "calc"):
        sys.exit("ERROR: mode must be gen | run | cal")
    xyz_list = cfg["xyz_list"]
    if not xyz_list:
        sys.exit("ERROR: no .xyz file(s) found in xyz_file of INPUT.cfg.")
    print(LOGO)
    if len(xyz_list) > 1:
        print(f"  Queued {len(xyz_list)} xyz files - running one by one:\n")
        for i, xyz in enumerate(xyz_list, 1):
            print(f"    [{i}] {xyz}")
        print()
    for xyz in xyz_list:
        run_one(xyz, cfg, args.verbose)

if __name__ == "__main__":
    main()
