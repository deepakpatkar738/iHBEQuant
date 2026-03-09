"""
iHBEQuant.py - Individual Hydrogen Bond Energy Quantifier
Author : Deepak Patkar
Usage  : python iHBEQuant.py [-c INPUT.cfg] [-v]
"""
import os, sys, math, configparser, argparse, re, subprocess
from datetime import datetime

G16    = "g16"           # Gaussian 16 executable
H2KCAL = 627.5094740631  # Hartree -> kcal/mol
COV_R  = {"H":0.23,"C":0.77,"N":0.75,"O":0.70,"F":0.71,
           "S":1.05,"P":1.06,"CL":0.99,"SI":1.11,"DEFAULT":0.90} # Covalent radii (Angstrom) for bond detection
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
║      iHBEQuant : Individual Hydrogen Bond Energy Quantifier | version 1.1             ║
║      Author    : Deepak Patkar                              | last update 2026        ║
╠═══════════════════════════════════════════════════════════════════════════════════════╣
║  Citations:                                                                           ║
║  1. Ahirwar MB; Gadre SR; Deshmukh MM J. Phys. Chem. A 2020, 124, 6699–6706.          ║
║  2. Patkar D; Ahirwar MB; Deshmukh MM ChemPhysChem 2022, 23, e202200476.              ║
║  3. Patkar D; Ahirwar MB; Deshmukh MM ChemPhysChem 2022, 23, e202200143.              ║
║  4. Patkar D; Ahirwar MB; Deshmukh MM New J. Chem. 2022, 46, 2368–2379.               ║
║  5. Patkar D; Ahirwar MB; Gadre SR; Deshmukh MM J. Phys. Chem. A 2021, 125, 8836–8845.║
╚═══════════════════════════════════════════════════════════════════════════════════════╝
"""
# ── config ────────────────────────────────────────────────────────────────────
def load_cfg(path):
    if not os.path.isfile(path): sys.exit(f"ERROR: config '{path}' not found.")
    # latin-1 safely reads any single-byte encoding (box chars, ¦, etc.)
    clean = []
    for l in open(path, encoding="latin-1"):
        l = re.sub(r'\s*#.*', '', l).rstrip()
        if not l.strip(): continue
        if l.strip().startswith("[") or "=" in l:
            if l.strip().startswith("[") and clean and clean[-1] != "": clean.append("")
            clean.append(l)
    c = configparser.ConfigParser(); c.read_string("\n".join(clean))
    MODE = {"gen":"generate","generate":"generate","run":"run","submit":"run","cal":"calc","calc":"calc"}
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
def load_xyz(path):                  # Read XYZ file 
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
                           "ha_dist":round(dha,3),"angle":round(ang,3)})
    return hbonds

# ── gaussian ──────────────────────────────────────────────────────────────────
def write_gjf(path,idx,atoms,cfg,title):
    extra = f" {cfg['keywords']}" if cfg.get("keywords") else ""
    with open(path,"w",encoding="utf-8") as f:
        f.write(f"%chk={os.path.splitext(path)[0]}.chk\n%nproc={cfg['nproc']}\n%mem={cfg['mem']}\n")
        f.write(f"#P {cfg['method']}/{cfg['basis']} SP nosymm{extra}\n\n{title}\n\n{cfg['charge']} {cfg['mult']}\n")
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

# ── fragments & energies ──────────────────────────────────────────────────────
def build_fragments(hb,n,monomers,atoms,cfg,out_dir,prefix):
    dm,am = hb["dm"],hb["am"]
    d_lbl,a_lbl = mono_label(dm,monomers,atoms),mono_label(am,monomers,atoms)
    tag,p = f"HB{n}",f"{out_dir}/{prefix}"
    all_i = [i for m in monomers for i in m]
    rm_d  = [i for mi,m in enumerate(monomers) for i in m if mi!=dm]
    rm_a  = [i for mi,m in enumerate(monomers) for i in m if mi!=am]
    rm_da = [i for mi,m in enumerate(monomers) for i in m if mi not in (dm,am)]
    files = {"cluster":f"{p}_cluster.gjf","rm_d":f"{p}_{tag}_rmD.gjf","rm_a":f"{p}_{tag}_rmA.gjf",
             "rm_da":f"{p}_{tag}_rmDA.gjf","dimer_da":f"{p}_{tag}_dimerDA.gjf",
             "dimer_d":f"{p}_{tag}_dimerD.gjf","dimer_a":f"{p}_{tag}_dimerA.gjf"}
    write_gjf(files["cluster"],  all_i,                     atoms,cfg,"Full cluster")
    write_gjf(files["rm_d"],     rm_d,                      atoms,cfg,f"{tag}_rmD | removed donor={d_lbl}")
    write_gjf(files["rm_a"],     rm_a,                      atoms,cfg,f"{tag}_rmA | removed acceptor={a_lbl}")
    write_gjf(files["rm_da"],    rm_da,                     atoms,cfg,f"{tag}_rmDA | removed D={d_lbl} A={a_lbl}")
    write_gjf(files["dimer_da"], monomers[dm]+monomers[am], atoms,cfg,f"{tag}_dimerDA | D={d_lbl} + A={a_lbl}")
    write_gjf(files["dimer_d"],  monomers[dm],              atoms,cfg,f"{tag}_dimerD | donor={d_lbl}")
    write_gjf(files["dimer_a"],  monomers[am],              atoms,cfg,f"{tag}_dimerA | acceptor={a_lbl}")
    return files

def compute_energies(files,mode):
    e = {k: get_energy(p,mode) for k,p in files.items()}
    res = {"hbe":None,"dimer":None,"coop":None}
    if mode == "generate": return res,e
    need = ["cluster","rm_d","rm_a","rm_da","dimer_da","dimer_d","dimer_a"]
    if any(e[k] is None for k in need): return res,e
    hbe_h   = abs(e["cluster"]-(e["rm_d"]+e["rm_a"]-e["rm_da"]))
    dimer_h = abs(e["dimer_d"]+e["dimer_a"]-e["dimer_da"])
    res["hbe"]   = round(hbe_h  *H2KCAL,3)
    res["dimer"] = round(dimer_h*H2KCAL,3)
    res["coop"]  = round((hbe_h-dimer_h)*H2KCAL,3)
    return res,e

def monomer_energies(atoms,monomers,cfg,out_dir,prefix,mode):
    mono_e = []
    for mi,m in enumerate(monomers):
        gjf = f"{out_dir}/{prefix}_M{mi+1}.gjf"
        write_gjf(gjf,m,atoms,cfg,f"Monomer M{mi+1} | {mono_label(mi,monomers,atoms)}")
        mono_e.append(get_energy(gjf,mode))
    return mono_e

def cross_validate(records):
    vals = [r["raw"].get("cluster") for r in records if r["raw"].get("cluster")]
    if len(set(vals)) > 1:
        print("  CROSS-VAL WARNING: cluster energy differs across HBs - check logs!")

# ── output ────────────────────────────────────────────────────────────────────
def fmt(v,dec=3): return f"{v:.{dec}f}" if v is not None else "pending"

def write_summary(path,records,atoms,monomers,mono_e,cluster_e,cfg,xyz,mode,run_time):
    n_hb,n_mono = len(records),len(monomers)
    T = "-"*95
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
    C1,C2 = 26,13
    lines = [LOGO,
             f"  File            : {xyz:<{C1}}  Mode           : {mode:<{C2}}  Run Time      : {run_time}",
             f"  Level of Theory : {lvl:<{C1}}  No. of monomers: {n_mono:<{C2}}  No. of HB     : {n_hb}",
             f"  Binding energy  : {fmt(abs(be_kcal) if be_kcal is not None else None):<{C1}}  Avg per HB     : {fmt(abs(avg_be) if avg_be is not None else None)}",
             f"  MTA_HBEs Sum    : {fmt(sum_hbe):<{C1}}",
             f"  Total MTA error : {fmt(total_err):<{C1}}  Per-HB error   : {fmt(per_hb_err)}",
             f"  # All energies in kcal/mol",""]

    W   = [5,18,10,10,11,15,15,15]
    SEP = "+"+"+" .join("-"*w for w in W)+"+"
    HDR = (f"|{'HB#':^{W[0]}}|{'D-H...A':^{W[1]}}|{'D-H(A)':^{W[2]}}"
           f"|{'H...A(A)':^{W[3]}}|{'Angle(deg)':^{W[4]}}"
           f"|{'MTA_HBE':^{W[5]}}|{'Dimer':^{W[6]}}|{'Coop':^{W[7]}}|")
    lines += [SEP,HDR,SEP]
    for r in records:
        hb,res = r["hb"],r["res"]
        label  = f"{lbl(hb['di'],atoms)}-{lbl(hb['hi'],atoms)}...{lbl(hb['ai'],atoms)}"
        lines.append(f"|{r['n']:^{W[0]}}|{label:^{W[1]}}|{fmt(hb['dh_dist']):^{W[2]}}"
                     f"|{fmt(hb['ha_dist']):^{W[3]}}|{fmt(hb['angle']):^{W[4]}}"
                     f"|{fmt(res['hbe']):^{W[5]}}|{fmt(res['dimer']):^{W[6]}}|{fmt(res['coop']):^{W[7]}}|")
    lines += [SEP,""]

    LBL = {"cluster":"Parent cluster","rm_d":"rmD","rm_a":"rmA","rm_da":"rmDA",
           "dimer_da":"dimerDA","dimer_d":"dimerD","dimer_a":"dimerA"}
    lines += ["  CALCULATION DETAILS","",
              f"  {xyz:<24} No of Atoms: {len(atoms):<8} No of Monomers: {n_mono:<8} HB detected: {n_hb}",""]
    for r in records:
        hb,res,raw = r["hb"],r["res"],r["raw"]
        tag   = f"HB{r['n']}"
        label = f"{lbl(hb['di'],atoms)}-{lbl(hb['hi'],atoms)}...{lbl(hb['ai'],atoms)}"
        lines += [f"  {tag}: {label}",T]
        for left,right in [("cluster","dimer_da"),("rm_a","dimer_a"),("rm_d","dimer_d"),("rm_da",None)]:
            ls = f"    {LBL[left]:<16}: {fmt(raw.get(left),6)}"
            rs = f"    {LBL[right]:<16}: {fmt(raw.get(right),6)}" if right else "                              # Note: energies in Hartree"
            lines.append(f"{ls:<45}{rs}")
        lines += [T, f"    MTA_{tag} : {fmt(res['hbe'])} kcal/mol   |   Dimer : {fmt(res['dimer'])} kcal/mol   |   Coop  : {fmt(res['coop'])} kcal/mol",""]

    text = "\n".join(lines)
    print(text)
    open(path,"w",encoding="utf-8").write(text)

def write_coords(path,records,atoms,monomers,prefix,cfg):
    FRAG_ORDER = ["dimer_da","dimer_d","dimer_a","rm_d","rm_a","rm_da"]
    FLBL = {"dimer_da":"dimerDA","dimer_d":"dimerD","dimer_a":"dimerA","rm_d":"rmD","rm_a":"rmA","rm_da":"rmDA"}
    def xyz_block(title,idx_list):
        lines = [str(len(idx_list)),title]
        for i in idx_list:
            a = atoms[i]; lines.append(f"{a['s']:<4}{a['x']:>14.8f}{a['y']:>14.8f}{a['z']:>14.8f}")
        return lines+[""]
    out = [f"# iHBQuant coordinates | {prefix}",f"# Level : {cfg['method']}/{cfg['basis']}",
           f"# Generated : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
           f"# Atoms : {len(atoms)}   Monomers : {len(monomers)}",""]
    out += xyz_block("Full cluster",[i for m in monomers for i in m])
    for r in records:
        hb  = r["hb"]; dm,am = hb["dm"],hb["am"]
        tag = f"HB{r['n']}: {lbl(hb['di'],atoms)}-{lbl(hb['hi'],atoms)}...{lbl(hb['ai'],atoms)}"
        frags = {"dimer_da":monomers[dm]+monomers[am],"dimer_d":monomers[dm],"dimer_a":monomers[am],
                 "rm_d":[i for mi,m in enumerate(monomers) for i in m if mi!=dm],
                 "rm_a":[i for mi,m in enumerate(monomers) for i in m if mi!=am],
                 "rm_da":[i for mi,m in enumerate(monomers) for i in m if mi not in (dm,am)]}
        out += [f"# {'='*70}",f"# {tag}",f"# {'='*70}",""]
        for key in FRAG_ORDER: out += xyz_block(f"{tag} | {FLBL[key]}",frags[key])
    open(path,"w",encoding="utf-8").write("\n".join(out))

# ── main ──────────────────────────────────────────────────────────────────────
def run_one(xyz,cfg,verbose):
    mode,run_time = cfg["mode"],datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    if not os.path.isfile(xyz): print(f"  ERROR: '{xyz}' not found - skipping.\n"); return
    prefix  = os.path.splitext(os.path.basename(xyz))[0]
    out_dir = prefix; os.makedirs(out_dir,exist_ok=True)
    atoms    = load_xyz(xyz); monomers = find_monomers(atoms); hbonds = find_hbonds(atoms,monomers,cfg)
    if not hbonds: print(f"  No HB found in {xyz}. Adjust thresholds in INPUT.cfg.\n"); return
    cluster_gjf = f"{out_dir}/{prefix}_cluster.gjf"
    write_gjf(cluster_gjf,[i for m in monomers for i in m],atoms,cfg,"Full cluster")
    cluster_e = get_energy(cluster_gjf,mode)
    mono_e    = monomer_energies(atoms,monomers,cfg,out_dir,prefix,mode)
    records   = []
    for n,hb in enumerate(hbonds,1):
        files,res_raw = build_fragments(hb,n,monomers,atoms,cfg,out_dir,prefix), None
        res,raw_e = compute_energies(files,mode)
        if verbose:
            for k,p in files.items(): print(f"    {os.path.basename(p)}")
        records.append({"n":n,"hb":hb,"res":res,"raw":raw_e})
    if mode in ("run","calc"): cross_validate(records)
    summary = f"{out_dir}/{prefix}_summary.txt"; coords = f"{out_dir}/{prefix}_coords.xyz"
    write_summary(summary,records,atoms,monomers,mono_e,cluster_e,cfg,xyz,mode,run_time)
    write_coords(coords,records,atoms,monomers,prefix,cfg)
    print(f"  Summary -> {summary}"); print(f"  Coords  -> {coords}\n")

def main():
    ap = argparse.ArgumentParser(description="iHBEQuant: Individual HB Energy Quantifier")
    ap.add_argument("-c","--config",  default="INPUT.cfg",help="config file [INPUT.cfg]")
    ap.add_argument("-v","--verbose", action="store_true", help="list gjf files written")
    args = ap.parse_args()
    cfg = load_cfg(args.config)
    if cfg["mode"] not in ("generate","run","calc"): sys.exit("ERROR: mode must be gen | run | cal")
    xyz_list = cfg["xyz_list"]
    if not xyz_list: sys.exit("ERROR: no .xyz file(s) found in xyz_file of INPUT.cfg.")
    print(LOGO)
    if len(xyz_list) > 1:
        print(f"  Queued {len(xyz_list)} xyz files - running one by one:\n")
        for i,xyz in enumerate(xyz_list,1): print(f"    [{i}] {xyz}")
        print()
    for xyz in xyz_list: run_one(xyz,cfg,args.verbose)

if __name__ == "__main__": main()
