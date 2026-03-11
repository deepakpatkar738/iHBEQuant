"""
iHBEQuant.py - Individual Hydrogen Bond Energy Quantifier
Author : Deepak Patkar
Version: 1.1
Usage  : python iHBEQuant.py [-c INPUT.cfg] [-v]

File-naming convention (example prefix = F1O2):
  Parent molecule  : F1O2/N1S2_M.xyz   F1O2/F1O2_M.gjf
  Monomers         : F1O2/N1S2_mono_M1.gjf  ...
  Primary frags    : F1O2/HB1_F1O2_frag_M1.gjf  HB1_F1O2_frag_M2.gjf
  Overlap frags    : F1O2/HB1_F1O2_frag_M1M2.gjf
  Dimer            : F1O2/HB1_F1O2_dimer_M1M2.gjf
  Summary          : F1O2/F1O2_summary.txt  (includeing coordinates)

Energies:
  MTA_HBn  = (E[frag_M1] + E[frag_M2] - E[frag_M1M2]) - E[parent_M]   (x 627.51)
  Dimer_HBn= (E[mono_M1] + E[mono_M2] - E[dimer_M1M2])                 (x 627.51)
  Coop_HBn = MTA_HBn - Dimer_HBn   #where n is HB index i.e., HB1, HB2, HB3...etc
"""
import os, sys, math, shutil, configparser, argparse, re, subprocess
from datetime import datetime

G16    = "g16"
H2KCAL = 627.5094740631
COV_R  = {"H":0.23,"C":0.77,"N":0.75,"O":0.70,"F":0.71,
           "S":1.05,"P":1.06,"CL":0.99,"SI":1.11,"DEFAULT":0.90}
HB_MTA = {"N","O","F","CL","S"}

LOGO = r"""
╔═══════════════════════════════════════════════════════════════════════════════════════╗
║                                                                                       ║
║        _   ____  ____  ______   ________    ___                              _        ║
║       (_) |_   ||   _||_   _ \ |_   __  | .'   `.                           / |_      ║
║       __    | |__| |    | |_) |  | |_ \_|/  .-.  \  __   _   ,--.   _ .--. `| |-'     ║
║      [  |   |  __  |    |  __'.  |  _| _ | |   | | [  | | | `'_\ : [ `.-. | | |       ║
║       | |  _| |  | |_  _| |__) |_| |__/ |\  `-'  \_ | \_/ |,// | |, | | | | | |,      ║
║      [___]|____||____||_______/|________| `.___.\__|'.__.'_/\'-;__/[___||__]\__/      ║
║                                                                                       ║
║      iHBEQuant : Individual Hydrogen Bond Energy Quantifier | version 1.2             ║
║      Author    : Deepak Patkar                              | last update 2026        ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║  Cite us:                                                                             ║
║  1. Patkar D; Ahirwar MB; Deshmukh MM ChemPhysChem 2022, 23, e202200476.              ║
║  2. Patkar D; Ahirwar MB; Deshmukh MM ChemPhysChem 2022, 23, e202200143.              ║
║  3. Patkar D; Ahirwar MB; Deshmukh MM New J. Chem. 2022, 46, 2368–2379.               ║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
"""

# ── config ────────────────────────────────────────────────────────────────────
def load_cfg(path):
    if not os.path.isfile(path): sys.exit(f"ERROR: config '{path}' not found.")
    clean = []
    for l in open(path, encoding="latin-1"):
        l = re.sub(r'\s*#.*', '', l).rstrip()
        if not l.strip(): continue
        if l.strip().startswith("[") or "=" in l:
            if l.strip().startswith("[") and clean and clean[-1] != "": clean.append("")
            clean.append(l)
    c = configparser.ConfigParser(); c.read_string("\n".join(clean))
    MODE = {"gen":"generate","generate":"generate","run":"run","submit":"run",
            "cal":"calc","calc":"calc"}
    xyz_raw = c.get("SYSTEM","xyz_file").strip()
    return {
        "xyz_list" : [x.strip() for x in xyz_raw.split() if x.strip().endswith(".xyz")],
        "nproc"    : c.getint("SYSTEM","nproc",              fallback=8),
        "mem"      : c.get("SYSTEM","mem",                   fallback="4GB").strip(),
        "method"   : c.get("THEORY","method",                fallback="B3LYP").strip(),
        "basis"    : c.get("THEORY","basis",                 fallback="6-311+G(d,p)").strip(),
        "keywords" : c.get("THEORY","keywords",              fallback="").strip(),
        "charge"   : c.getint("THEORY","charge",             fallback=0),
        "mult"     : c.getint("THEORY","multiplicity",       fallback=1),
        "d_min"    : c.getfloat("THRESHOLDS","hb_min_dist",  fallback=1.2),
        "d_max"    : c.getfloat("THRESHOLDS","hb_max_dist",  fallback=2.5),
        "a_min"    : c.getfloat("THRESHOLDS","hb_min_angle", fallback=100.0),
        "a_max"    : c.getfloat("THRESHOLDS","hb_max_angle", fallback=180.0),
        "nb_th"    : c.getfloat("THRESHOLDS","neighbor_th",  fallback=2.5),
        "mode"     : MODE.get(c.get("RUN","mode",fallback="gen").strip().lower(),"generate"),
    }

# ── geometry ──────────────────────────────────────────────────────────────────
def load_xyz(path):
    raw = open(path).readlines()
    def _parse(lines):
        out = []
        for l in lines:
            p = l.split()
            if len(p) < 4 or p[0].lstrip("+-").replace(".","").isdigit(): continue
            try: x,y,z = float(p[1]),float(p[2]),float(p[3])
            except ValueError: continue
            out.append({"s":p[0],"x":x,"y":y,"z":z,
                        "gjf":f"{p[0]:<4}{x:>14.8f}{y:>14.8f}{z:>14.8f}\n"})
        return out
    try:
        n = int(raw[0].strip()); return _parse(raw[2:2+n])
    except (ValueError, IndexError): return _parse(raw)

def dist(a,b):    return math.sqrt((a["x"]-b["x"])**2+(a["y"]-b["y"])**2+(a["z"]-b["z"])**2)
def lbl(i,atoms): return f"{atoms[i]['s']}{i+1}"
def mono_label(mi,monomers,atoms): return "-".join(lbl(a,atoms) for a in monomers[mi])

def bonded(a,b,tol=0.5):
    r = COV_R.get(a["s"].upper(),COV_R["DEFAULT"]) + COV_R.get(b["s"].upper(),COV_R["DEFAULT"])
    return abs(dist(a,b)-r) <= tol

def hb_angle(donor,h,acc):
    def v(o,t): return [t["x"]-o["x"],t["y"]-o["y"],t["z"]-o["z"]]
    v1,v2 = v(h,donor),v(h,acc)
    dot = sum(i*j for i,j in zip(v1,v2))
    nm  = lambda u: math.sqrt(sum(i*i for i in u))
    return math.degrees(math.acos(max(-1.0,min(1.0,dot/(nm(v1)*nm(v2)+1e-12)))))

def find_monomers(atoms):
    visited,monomers = set(),[]
    def fill(start):
        stack,m = [start],[]
        while stack:
            i = stack.pop()
            if i in visited: continue
            visited.add(i); m.append(i)
            for j in range(len(atoms)):
                if j not in visited and bonded(atoms[i],atoms[j]): stack.append(j)
        return m
    for i in range(len(atoms)):
        if i not in visited: monomers.append(fill(i))
    return monomers

def find_hbonds(atoms,monomers,cfg):
    mono_of = {a:mi for mi,m in enumerate(monomers) for a in m}
    hbonds  = []
    for hi in range(len(atoms)):
        if atoms[hi]["s"].upper() != "H": continue
        donor = next((d for d in range(len(atoms))
                      if d!=hi and atoms[d]["s"].upper()!="H" and bonded(atoms[hi],atoms[d])),None)
        if donor is None: continue
        for ai in range(len(atoms)):
            if ai in (hi,donor) or atoms[ai]["s"].upper() not in HB_MTA: continue
            if mono_of[ai] == mono_of[hi]: continue
            dha = dist(atoms[hi],atoms[ai])
            if not (cfg["d_min"] < dha < cfg["d_max"]): continue
            ang = hb_angle(atoms[donor],atoms[hi],atoms[ai])
            if not (cfg["a_min"] <= ang <= cfg["a_max"]): continue
            hbonds.append({"di":donor,"hi":hi,"ai":ai,"dm":mono_of[donor],"am":mono_of[ai],
                           "dh_dist":round(dist(atoms[donor],atoms[hi]),3),
                           "ha_dist":round(dha,3),"angle":round(ang,3),
                           "da_dist":round(dist(atoms[donor],atoms[ai]),3)})
    return hbonds

# ── gaussian ──────────────────────────────────────────────────────────────────
def route_line(cfg):
    """Build the Gaussian route line, appending keywords if provided."""
    kw = f" {cfg['keywords'].strip()}" if cfg.get("keywords","").strip() else ""
    return f"#P {cfg['method']}/{cfg['basis']} SP nosymm{kw}"

def write_gjf(path,idx,atoms,cfg,title):
    with open(path,"w",encoding="utf-8") as f:
        f.write(f"%chk={os.path.splitext(path)[0]}.chk\n"
                f"%nproc={cfg['nproc']}\n%mem={cfg['mem']}\n")
        f.write(f"{route_line(cfg)}\n\n{title}\n\n{cfg['charge']} {cfg['mult']}\n")
        for i in idx: f.write(atoms[i]["gjf"])
        f.write("\n\n")

def log_path(gjf):   return gjf.replace(".gjf",".log")
def normal_end(log): return os.path.isfile(log) and any("Normal termination" in l for l in open(log))

def run_g16(gjf):
    if normal_end(log_path(gjf)): return
    print(f"  Running g16: {os.path.basename(gjf)}")
    ret = subprocess.run([G16,gjf],capture_output=False)
    if ret.returncode != 0: print(f"  WARNING: g16 non-zero exit for {gjf}")

def parse_energy(log):
    e = None
    for l in open(log):
        if "SCF Done"  in l: e = float(l.split()[4])
        if "EUMP2"     in l: e = float(l.split()[5].replace("D","E"))
        if "E(Corr)="  in l: e = float(l.split()[4])
        if "CCSD(T)= " in l: e = float(l.split()[1].replace("D","E"))
    return e

def get_energy(gjf,mode):
    if mode == "generate": return None
    if mode == "run": run_g16(gjf)
    lp = log_path(gjf)
    if normal_end(lp): return parse_energy(lp)
    print(f"  WARNING: abnormal/missing log - {gjf}"); return None

# ── naming helpers ─────────────────────────────────────────────────────────────
def parent_gjf(out_dir,prefix):
    """N1S2/N1S2_M.gjf"""
    return f"{out_dir}/{prefix}_M.gjf"

def mono_gjf(out_dir,prefix,mi):
    """N1S2/N1S2_mono_M1.gjf  (mi is 0-based)"""
    return f"{out_dir}/{prefix}_mono_M{mi+1}.gjf"

def frag_gjf(out_dir,prefix,n,dm,am,kind):
    """
    kind='donor'   -> HB1_N1S2_frag_M{dm+1}.gjf    (primary fragment: donor side)
    kind='acceptor'-> HB1_N1S2_frag_M{am+1}.gjf    (primary fragment: acceptor side)
    kind='overlap' -> HB1_N1S2_frag_M{dm+1}M{am+1}.gjf  (overlap fragment)
    """
    tag = f"HB{n}_{prefix}"
    if   kind == "donor":    return f"{out_dir}/{tag}_frag_M{dm+1}.gjf"
    elif kind == "acceptor": return f"{out_dir}/{tag}_frag_M{am+1}.gjf"
    elif kind == "overlap":  return f"{out_dir}/{tag}_frag_M{dm+1}M{am+1}.gjf"

def dimer_gjf(out_dir,prefix,n,dm,am):
    """HB1_N1S2_dimer_M1M2.gjf"""
    return f"{out_dir}/HB{n}_{prefix}_dimer_M{dm+1}M{am+1}.gjf"

# ── fragments & energies ──────────────────────────────────────────────────────
def build_fragments(hb,n,monomers,atoms,cfg,out_dir,prefix):
    dm,am = hb["dm"],hb["am"]
    d_lbl,a_lbl = mono_label(dm,monomers,atoms),mono_label(am,monomers,atoms)
    all_i  = [i for m in monomers for i in m]
    # primary fragments: cluster minus one monomer
    frag_d_idx = [i for mi,m in enumerate(monomers) for i in m if mi!=dm]   # remove donor
    frag_a_idx = [i for mi,m in enumerate(monomers) for i in m if mi!=am]   # remove acceptor
    # overlap fragment: cluster minus both donor & acceptor
    frag_da_idx= [i for mi,m in enumerate(monomers) for i in m if mi not in (dm,am)]
    # dimer: only the donor-acceptor pair
    dimer_idx  = monomers[dm]+monomers[am]

    p_gjf  = parent_gjf(out_dir,prefix)
    fd_gjf = frag_gjf(out_dir,prefix,n,dm,am,"donor")
    fa_gjf = frag_gjf(out_dir,prefix,n,dm,am,"acceptor")
    fo_gjf = frag_gjf(out_dir,prefix,n,dm,am,"overlap")
    di_gjf = dimer_gjf(out_dir,prefix,n,dm,am)

    # parent molecule written once (reused across HBs)
    if not os.path.isfile(p_gjf):
        write_gjf(p_gjf, all_i, atoms, cfg,
                  f"{prefix} parent molecule | {len(atoms)} atoms {len(monomers)} monomers")

    write_gjf(fd_gjf, frag_d_idx, atoms, cfg,
              f"HB{n} primary frag (donor side) | removed acceptor M{am+1}={a_lbl}")
    write_gjf(fa_gjf, frag_a_idx, atoms, cfg,
              f"HB{n} primary frag (acceptor side) | removed donor M{dm+1}={d_lbl}")
    write_gjf(fo_gjf, frag_da_idx, atoms, cfg,
              f"HB{n} overlap frag | removed D={d_lbl} A={a_lbl}")
    write_gjf(di_gjf, dimer_idx, atoms, cfg,
              f"HB{n} dimer M{dm+1}M{am+1} | D={d_lbl} + A={a_lbl}")

    return {"parent":p_gjf,"frag_d":fd_gjf,"frag_a":fa_gjf,
            "frag_da":fo_gjf,"dimer":di_gjf,
            # monomer keys filled by caller
            "mono_d":mono_gjf(out_dir,prefix,dm),
            "mono_a":mono_gjf(out_dir,prefix,am)}

def compute_energies(files,mode):
    """
    MTA_HBn  = ( E[frag_d] + E[frag_a] - E[frag_da] ) - E[parent]   × 627.51
    Dimer_HBn= ( E[mono_d] + E[mono_a] - E[dimer]   )               × 627.51
    Coop_HBn = MTA_HBn - Dimer_HBn
    """
    e = {k: get_energy(p,mode) for k,p in files.items()}
    res = {"hbe":None,"dimer":None,"coop":None}
    if mode == "generate": return res,e
    need = ["parent","frag_d","frag_a","frag_da","dimer","mono_d","mono_a"]
    if any(e.get(k) is None for k in need): return res,e
    mta_h   = (e["frag_d"] + e["frag_a"] - e["frag_da"]) - e["parent"]
    dimer_h = (e["mono_d"] + e["mono_a"] - e["dimer"])
    res["hbe"]   = round(mta_h  * H2KCAL, 3)
    res["dimer"] = round(dimer_h* H2KCAL, 3)
    res["coop"]  = round((mta_h - dimer_h) * H2KCAL, 3)
    return res,e

def monomer_energies(atoms,monomers,cfg,out_dir,prefix,mode):
    mono_e = []
    for mi,m in enumerate(monomers):
        gjf = mono_gjf(out_dir,prefix,mi)
        write_gjf(gjf,m,atoms,cfg,
                  f"{prefix} monomer M{mi+1} | {mono_label(mi,monomers,atoms)}")
        mono_e.append(get_energy(gjf,mode))
    return mono_e

def cross_validate(records):
    vals = [r["raw"].get("parent") for r in records if r["raw"].get("parent")]
    if len(set(vals)) > 1:
        print("  CROSS-VAL WARNING: parent energy differs across HBs - check logs!")

# ── output ────────────────────────────────────────────────────────────────────
def fmt(v,dec=3): return f"{v:.{dec}f}" if v is not None else "pending"

def xyz_block_str(title,idx,atoms):
    """Return an XYZ-formatted block as a string."""
    lines = [str(len(idx)), title]
    for i in idx:
        a = atoms[i]
        lines.append(f"{a['s']:<4}{a['x']:>14.8f}{a['y']:>14.8f}{a['z']:>14.8f}")
    lines.append("")
    return "\n".join(lines)

def write_summary(path,records,atoms,monomers,mono_e,cluster_e,cfg,xyz,mode,run_time,prefix,out_dir):
    n_hb,n_mono = len(records),len(monomers)
    T = "-"*95

    # binding energy block
    if None not in [cluster_e]+mono_e and mode in ("run","calc"):
        be_kcal = (cluster_e-sum(mono_e))*H2KCAL; avg_be = be_kcal/n_hb if n_hb else None
    else: be_kcal = avg_be = None
    hbe_vals  = [r["res"]["hbe"] for r in records if r["res"]["hbe"] is not None]
    sum_hbe   = round(sum(hbe_vals),3) if hbe_vals else None
    if None not in (be_kcal,sum_hbe):
        total_err  = round(abs(abs(be_kcal)-sum_hbe),3)
        per_hb_err = round(total_err/n_hb,3) if n_hb else None
    else: total_err = per_hb_err = None

    lvl = f"{cfg['method']}/{cfg['basis']}"
    kw_str = cfg.get("keywords","").strip()
    if kw_str: lvl = f"{lvl}  {kw_str}"
    avg_mta = round(sum_hbe/n_hb,3) if (sum_hbe is not None and n_hb) else None

    def hrow(ll, lv, rl="", rv=""):
        cell  = f"  {ll:<18}: {str(lv)}"
        left  = f"{cell:<56}"
        if not rl:   right = ""
        elif not rv: right = f"  {rl}"          # comment line — no colon
        else:        right = f"  {rl:<20}: {rv}"
        return left + right

    lines = [LOGO,
        hrow("Run Time",        run_time,        "No. of monomers", n_mono),
        hrow("File",            xyz,             "No. of HB",       n_hb),
        hrow("Mode",            mode,            "Binding energy",  fmt(abs(be_kcal) if be_kcal is not None else None)),
        hrow("Level of Theory", lvl,             "Avg per HB",      fmt(abs(avg_be)  if avg_be  is not None else None)),
        hrow("hb_min_dist",     cfg["d_min"],    "MTA_HBEs Sum",    fmt(sum_hbe)),
        hrow("hb_max_dist",     cfg["d_max"],    "Avg per MTA_HB",  fmt(avg_mta)),
        hrow("hb_min_angle",    cfg["a_min"],    "Total MTA error", fmt(total_err)),
        hrow("hb_max_angle",    cfg["a_max"],    "Per-HB error",    fmt(per_hb_err)),
        hrow("neighbor_th",     cfg["nb_th"],    "# All energies in kcal/mol", ""),
        ""]

    # HB table
    W   = [5,18,10,10,10,11,15,15,15]
    SEP = "+"+"+" .join("-"*w for w in W)+"+"
    HDR = (f"|{'HB#':^{W[0]}}|{'D-H...A':^{W[1]}}|{'D-H(A)':^{W[2]}}"
           f"|{'H...A(A)':^{W[3]}}|{'A...D(A)':^{W[4]}}|{'Angle(deg)':^{W[5]}}"
           f"|{'MTA_HBE':^{W[6]}}|{'Dimer':^{W[7]}}|{'Coop':^{W[8]}}|")
    lines += [SEP,HDR,SEP]
    for r in records:
        hb,res = r["hb"],r["res"]
        label  = f"{lbl(hb['di'],atoms)}-{lbl(hb['hi'],atoms)}...{lbl(hb['ai'],atoms)}"
        lines.append(f"|{r['n']:^{W[0]}}|{label:^{W[1]}}|{fmt(hb['dh_dist']):^{W[2]}}"
                     f"|{fmt(hb['ha_dist']):^{W[3]}}|{fmt(hb['da_dist']):^{W[4]}}|{fmt(hb['angle']):^{W[5]}}"
                     f"|{fmt(res['hbe']):^{W[6]}}|{fmt(res['dimer']):^{W[7]}}|{fmt(res['coop']):^{W[8]}}|")
    lines += [SEP,""]

    # calculation details
    lines += ["  CALCULATION DETAILS","",
              f"  {xyz:<24} No of Atoms: {len(atoms):<8} No of Monomers: {n_mono:<8} HB detected: {n_hb}",""]

    for r in records:
        hb,res,raw = r["hb"],r["res"],r["raw"]
        dm,am = hb["dm"],hb["am"]
        tag   = f"HB{r['n']}"
        label = f"{lbl(hb['di'],atoms)}-{lbl(hb['hi'],atoms)}...{lbl(hb['ai'],atoms)}"
        lines += [f"  {tag}: {label}",T]
        # two-column energy table: label + value only, no filenames
        def erow(ll, lv, rl=None, rv=None):
            left  = f"    {ll:<12}  {fmt(lv,6)}"
            right = f"    {rl:<12}  {fmt(rv,6)}" if rl else "                              # Hartree"
            lines.append(f"{left:<42}{right}")
        erow("parent_M",   raw.get("parent"),  "dimer_M1M2", raw.get("dimer"))
        erow(f"frag_M{dm+1}",     raw.get("frag_d"),  f"mono_M{dm+1}",    raw.get("mono_d"))
        erow(f"frag_M{am+1}",     raw.get("frag_a"),  f"mono_M{am+1}",    raw.get("mono_a"))
        erow(f"frag_M{dm+1}M{am+1}", raw.get("frag_da"))

        lines += [T,
                  f"    MTA_{tag}  : {fmt(res['hbe'])} kcal/mol   |   Dimer : {fmt(res['dimer'])} kcal/mol   |   Coop  : {fmt(res['coop'])} kcal/mol",
                  ""]

    # ── Coordinates section (inline, no separate file) ────────────────────────
    lines += ["","  COORDINATES","  " + "="*91,""]
    all_i = [i for m in monomers for i in m]
    lines.append(f"  # Parent molecule : {prefix}_M  ({len(all_i)} atoms, {n_mono} monomers)")
    lines.append(xyz_block_str(f"{prefix} parent molecule",all_i,atoms))

    for mi,m in enumerate(monomers):
        lines.append(f"  # Monomer M{mi+1} : {mono_label(mi,monomers,atoms)}")
        lines.append(xyz_block_str(f"{prefix}_mono_M{mi+1}",m,atoms))

    FRAG_ORDER = [
        ("frag_donor",    "primary frag (donor side)"),
        ("frag_acceptor", "primary frag (acceptor side)"),
        ("frag_overlap",  "overlap frag"),
        ("dimer",         "dimer"),
    ]
    for r in records:
        hb  = r["hb"]; dm,am = hb["dm"],hb["am"]
        tag = f"HB{r['n']}: {lbl(hb['di'],atoms)}-{lbl(hb['hi'],atoms)}...{lbl(hb['ai'],atoms)}"
        lines += [f"  # {'='*70}",f"  # {tag}",f"  # {'='*70}",""]
        frag_map = {
            "frag_donor":    [i for mi,m in enumerate(monomers) for i in m if mi!=dm],
            "frag_acceptor": [i for mi,m in enumerate(monomers) for i in m if mi!=am],
            "frag_overlap":  [i for mi,m in enumerate(monomers) for i in m if mi not in (dm,am)],
            "dimer":         monomers[dm]+monomers[am],
        }
        for key,desc in FRAG_ORDER:
            fname = {
                "frag_donor":    os.path.basename(frag_gjf(out_dir,prefix,r['n'],dm,am,"donor")),
                "frag_acceptor": os.path.basename(frag_gjf(out_dir,prefix,r['n'],dm,am,"acceptor")),
                "frag_overlap":  os.path.basename(frag_gjf(out_dir,prefix,r['n'],dm,am,"overlap")),
                "dimer":         os.path.basename(dimer_gjf(out_dir,prefix,r['n'],dm,am)),
            }[key]
            lines.append(f"  # {tag} | {desc} [{fname}]")
            lines.append(xyz_block_str(f"{tag} | {desc}",frag_map[key],atoms))

    text = "\n".join(lines)
    print(text)
    open(path,"w",encoding="utf-8").write(text)

# ── main ──────────────────────────────────────────────────────────────────────
def run_one(xyz,cfg,verbose):
    mode     = cfg["mode"]
    run_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not os.path.isfile(xyz):
        print(f"  ERROR: '{xyz}' not found - skipping.\n"); return

    prefix  = os.path.splitext(os.path.basename(xyz))[0]
    out_dir = prefix
    os.makedirs(out_dir, exist_ok=True)

    # copy xyz into directory as <prefix>_M.xyz (parent molecule)
    parent_xyz = f"{out_dir}/{prefix}_M.xyz"
    shutil.copy2(xyz, parent_xyz)
    print(f"  Parent molecule  : {parent_xyz}")

    atoms    = load_xyz(xyz)
    monomers = find_monomers(atoms)
    hbonds   = find_hbonds(atoms,monomers,cfg)

    print(f"  Monomers detected: {len(monomers)}")
    print(f"  H-bonds detected : {len(hbonds)}")

    if not hbonds:
        print(f"  No HB found in {xyz}. Adjust thresholds in INPUT.cfg.\n"); return

    # write parent molecule gjf
    p_gjf = parent_gjf(out_dir,prefix)
    write_gjf(p_gjf, [i for m in monomers for i in m], atoms, cfg,
              f"{prefix} parent molecule | {len(atoms)} atoms {len(monomers)} monomers")
    cluster_e = get_energy(p_gjf, mode)

    # write monomers
    mono_e = monomer_energies(atoms,monomers,cfg,out_dir,prefix,mode)

    # process each HB
    records = []
    for n,hb in enumerate(hbonds,1):
        files   = build_fragments(hb,n,monomers,atoms,cfg,out_dir,prefix)
        res,raw_e = compute_energies(files,mode)
        if verbose:
            print(f"    HB{n} files:")
            for k,p in files.items(): print(f"      {k:<12}: {os.path.basename(p)}")
        records.append({"n":n,"hb":hb,"res":res,"raw":raw_e})

    if mode in ("run","calc"): cross_validate(records)

    summary = f"{out_dir}/{prefix}_summary.txt"
    write_summary(summary,records,atoms,monomers,mono_e,cluster_e,cfg,xyz,mode,run_time,prefix,out_dir)
    print(f"\n  Summary -> {summary}  (coordinates included)\n")

def main():
    ap = argparse.ArgumentParser(description="iHBEQuant: Individual HB Energy Quantifier")
    ap.add_argument("-c","--config",  default="INPUT.cfg", help="config file [INPUT.cfg]")
    ap.add_argument("-v","--verbose", action="store_true",  help="list gjf files written")
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    if cfg["mode"] not in ("generate","run","calc"):
        sys.exit("ERROR: mode must be gen | run | cal")
    xyz_list = cfg["xyz_list"]
    if not xyz_list: sys.exit("ERROR: no .xyz file(s) found in xyz_file of INPUT.cfg.")
    print(LOGO)
    if len(xyz_list) > 1:
        print(f"  Queued {len(xyz_list)} xyz files - running one by one:\n")
        for i,xyz in enumerate(xyz_list,1): print(f"    [{i}] {xyz}")
        print()
    for xyz in xyz_list: run_one(xyz,cfg,args.verbose)

if __name__ == "__main__": main()
